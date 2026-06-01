#!/usr/bin/env python
# Generated with assist from Claude code (Opus 4.8)
"""Plan (and optionally run) per-tissue post-processing of DSA-aligned PacBio CRAMs.

Driven by the decode table produced by make-manifest.py. Rows are grouped by a
key (default: the tissue code parsed from the sample name), each sample is mapped
to its DSA-aligned result CRAM, and per group the tool:

  1. merges the group's CRAMs into one (single-member groups are copied),
  2. runs mosdepth coverage on the merged CRAM (optional, --no-mosdepth),
  3. extracts CpG (mCG) scores with pb-CpG-tools (optional, --no-mcg).

By default this only writes a plan and prints the commands; pass --run to execute
them locally (heavy -- prefer the generated SLURM scripts via --slurm).
"""

import logging
import shutil
import subprocess
import sys
from pathlib import Path

import defopt
import pandas as pd

RESULT_GLOB = "{sample}-*_DSA.aligned.sorted.cram"

PLAN_COLUMNS = ["group", "action", "status", "n_inputs", "output", "inputs"]

# Default pb-CpG-tools binary (aligned_bam_to_cpg_scores).
DEFAULT_BAMTOCPG = (
    "/mmfs1/gscratch/stergachislab/mhsohny/Tools/"
    "pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores"
)
MOSDEPTH_MAX_THREADS = 4


def group_value(row, key):
    """Resolve one grouping key for a decode row. The special key 'tissue_code'
    is parsed from field 1 of the sample name; anything else is a column name."""
    if key == "tissue_code":
        return row["sample"].split("-")[1]
    if key not in row:
        raise ValueError(f"--group-by key {key!r} is not a decode column")
    return str(row[key])


def find_result_cram(results_dir, sample):
    """Locate the single DSA-aligned result CRAM for a sample, or None if absent."""
    matches = [
        m
        for m in results_dir.glob(RESULT_GLOB.format(sample=sample))
        if m.suffix == ".cram"
    ]
    if len(matches) > 1:
        raise ValueError(
            f"Sample {sample!r}: {len(matches)} result CRAMs match "
            f"{RESULT_GLOB.format(sample=sample)!r}: {sorted(m.name for m in matches)}"
        )
    return matches[0] if matches else None


def aligner_token(cram_name, sample):
    """Extract the '<aligner>_<version>' token from a result CRAM filename."""
    # {sample}-{aligner}_DSA.aligned.sorted.cram
    rest = cram_name[len(sample) + 1 :]
    return rest.split("_DSA.aligned.sorted.cram")[0]


def build_plan(decode_tsv, results_dir, output_dir, group_keys, skip_merge=False):
    """Group decode rows and resolve each group's inputs/output. Returns a list of
    dicts describing the planned merges/copies."""
    df = pd.read_csv(decode_tsv, sep="\t", dtype=str)

    groups = {}  # label -> list of row dicts
    for row in df.to_dict("records"):
        label = "_".join(group_value(row, k) for k in group_keys)
        groups.setdefault(label, []).append(row)

    plan = []
    for label, rows in sorted(groups.items()):
        donors = {r["donor"] for r in rows}
        dsas = {r["dsa"] for r in rows}
        if len(donors) > 1:
            raise ValueError(f"Group {label!r} spans multiple donors: {sorted(donors)}")
        if len(dsas) > 1:
            raise ValueError(f"Group {label!r} spans multiple DSAs: {sorted(dsas)}")
        donor = rows[0]["donor"]
        dsa = rows[0]["dsa"]

        inputs = []
        missing = []
        tokens = set()
        for r in rows:
            cram = find_result_cram(results_dir, r["sample"])
            if cram is None:
                missing.append(r["sample"])
            else:
                inputs.append(cram)
                tokens.add(aligner_token(cram.name, r["sample"]))

        status = "ready" if not missing else f"missing:{len(missing)}"
        if len(tokens) > 1:
            token = "mixed"
            status = "error:mixed-aligner"
        elif len(tokens) == 1:
            token = next(iter(tokens))
        else:
            token = "NA"  # no inputs found yet (all missing)
        # Name the merged CRAM with the tissue alias (e.g. "liver") when the group
        # shares one; fall back to the grouping label otherwise.
        aliases = {r["tissue_alias"] for r in rows}
        name = aliases.pop() if len(aliases) == 1 else label
        output = output_dir / f"{donor}-{name}-{token}_DSA.aligned.sorted.cram"
        action = "copy" if len(rows) == 1 else "merge"

        # With --skip-merge we run mosdepth/mCG on an already-merged CRAM, so
        # readiness depends on that output existing, not on the per-cell inputs.
        if skip_merge:
            status = "ready" if output.exists() else "missing-merged"

        plan.append(
            {
                "group": label,
                "action": action,
                "status": status,
                "n_inputs": len(rows),
                "output": str(output),
                "inputs": [str(p) for p in inputs],
                "dsa": dsa,
                "missing": missing,
                "donor": donor,
                "name": name,
                "samples": [r["sample"] for r in rows],
            }
        )
    return plan


def merge_command(entry, threads):
    """samtools merge command for a multi-member group.

    Inputs are embed_ref=1 CRAMs, so no external --reference is needed.
    """
    return [
        "samtools",
        "merge",
        "-@",
        str(threads),
        "-O",
        "CRAM",
        "--write-index",
        "--output-fmt-option",
        "embed_ref=1",
        "--output-fmt-option",
        "store_md=1",
        "--output-fmt-option",
        "store_nm=1",
        "-o",
        entry["output"],
        *entry["inputs"],
    ]


def mosdepth_commands(merged_cram, dsa, threads):
    """mosdepth coverage on the merged CRAM (--no-per-base; capped at 4 threads)."""
    parent = Path(merged_cram).parent
    depth_dir = parent / "Depth"
    prefix = depth_dir / Path(merged_cram).name.replace(".cram", "")
    return [
        ["mkdir", "-p", str(depth_dir)],
        [
            "mosdepth",
            "--threads",
            str(min(MOSDEPTH_MAX_THREADS, threads)),
            "--fasta",
            dsa,
            "--no-per-base",
            "--mapq",
            "0",
            str(prefix),
            merged_cram,
        ],
    ]


def mcg_commands(merged_cram, dsa, threads, bamtocpg):
    """pb-CpG-tools mCG extraction on the merged CRAM, haplotype-split via HP."""
    parent = Path(merged_cram).parent
    mcg_dir = parent / "mCG"
    prefix = mcg_dir / Path(merged_cram).name.replace(".cram", "")
    return [
        ["mkdir", "-p", str(mcg_dir)],
        [
            bamtocpg,
            "--bam",
            merged_cram,
            "--ref",
            dsa,
            "--output-prefix",
            str(prefix),
            "--pileup-mode",
            "model",
            "--modsites-mode",
            "denovo",
            "--min-coverage",
            "4",
            "--min-mapq",
            "1",
            "--hap-tag",
            "HP",
            "--threads",
            str(threads),
        ],
    ]


def group_commands(entry, threads, *, mosdepth, mcg, bamtocpg, skip_merge=False):
    """Full command sequence for one group: merge/copy, then optional mosdepth and
    mCG extraction on the resulting merged CRAM. With skip_merge, the merge/copy
    line is kept but commented out (the merged CRAM is assumed to exist already).

    Commands are arg lists; a leading "#" token marks a command as commented (not
    executed) when rendered to the script.
    """
    if entry["action"] == "copy":
        merge_cmds = [
            ["cp", entry["inputs"][0], entry["output"]],
            ["cp", entry["inputs"][0] + ".crai", entry["output"] + ".crai"],
        ]
    else:
        merge_cmds = [merge_command(entry, threads)]
    if skip_merge:
        merge_cmds = [["#", *c] for c in merge_cmds]

    cmds = list(merge_cmds)
    if mosdepth:
        cmds += mosdepth_commands(entry["output"], entry["dsa"], threads)
    if mcg:
        cmds += mcg_commands(entry["output"], entry["dsa"], threads, bamtocpg)
    return cmds


def render_slurm_script(
    job_name, body, *, threads, account, partition, time, mem_gb, mail
):
    """Render a standalone SLURM batch script (header + logging + set -euo
    pipefail) wrapping `body`, following the lab template."""
    return f"""#!/bin/bash
# generated by merge-by-tissue-pacbio-n-extract-mcg.py

#SBATCH --job-name={job_name}
#SBATCH --account={account}
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --cpus-per-task={threads}
#SBATCH --time={time}
#SBATCH --mem={mem_gb}G
#SBATCH -o logs/%x.%N.%j.slurm.stdout.log
#SBATCH -e logs/%x.%N.%j.slurm.stderr.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user={mail}
#SBATCH --export=ALL

set -euo pipefail
mkdir -p logs

echo "$(date '+%m-%d-%Y %H:%M:%S') JOBID: ${{SLURM_JOB_ID}} CMD: sbatch $(basename "$0")" >> logs/job_submission_commandline.log

{body}
"""


def main(
    decode_tsv: Path,
    *,
    results_dir: Path = Path("results"),
    output_dir: Path = Path("results"),
    group_by: str = "tissue_code",
    threads: int = 39,
    mosdepth: bool = True,
    mcg: bool = True,
    skip_merge: bool = False,
    bamtocpg: str = DEFAULT_BAMTOCPG,
    plan: Path = Path("merge_plan.tsv"),
    run: bool = False,
    slurm: bool = False,
    unified: bool = False,
    slurm_dir: Path = Path("slurm_merge"),
    account: str = "stergachislab",
    partition: str = "cpu-g2",
    slurm_time: str = "24:00:00",
    mem_gb: int = 32,
    mail: str = "mhsohny@uw.edu",
    verbose: int = 0,
):
    """
    Per-tissue post-processing of DSA-aligned PacBio CRAMs from a decode table:
    merge, then mosdepth coverage and pb-CpG-tools mCG extraction.

    :param decode_tsv: Decode table from make-manifest.py (tab-separated, with header)
    :param results_dir: Directory holding the {sample}-..._DSA.aligned.sorted.cram files
    :param output_dir: Directory to write merged CRAMs (and Depth/, mCG/) into
    :param group_by: Comma-separated grouping keys (decode columns or 'tissue_code')
    :param threads: Threads for merge + pb-CpG and --cpus-per-task (mosdepth capped at 4)
    :param mosdepth: Run mosdepth coverage on each merged CRAM
    :param mcg: Run pb-CpG-tools mCG extraction on each merged CRAM
    :param skip_merge: Skip the merge/copy step and run only mosdepth/mCG on
        already-merged CRAMs in --output-dir
    :param bamtocpg: Path to the pb-CpG-tools aligned_bam_to_cpg_scores binary
    :param plan: Output path for the merge-plan table
    :param run: Execute the planned commands locally (heavy; prefer the cluster)
    :param slurm: Emit SLURM batch scripts instead of printing raw commands
    :param unified: With --slurm, emit ONE script running all groups sequentially
        (instead of one parallel job per group)
    :param slurm_dir: Directory to write the SLURM script(s) into
    :param account: SLURM account for the batch scripts
    :param partition: SLURM partition for the batch scripts
    :param slurm_time: SLURM walltime (HH:MM:SS) for the batch scripts
    :param mem_gb: Memory (GB) requested per SLURM job
    :param mail: Email for SLURM --mail-user
    :param verbose: Set the logging level of the function
    """
    logging.basicConfig(
        format="[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    )
    logging.getLogger().setLevel(10 * (3 - verbose))

    group_keys = [k.strip() for k in group_by.split(",") if k.strip()]
    entries = build_plan(
        decode_tsv, results_dir, output_dir, group_keys, skip_merge=skip_merge
    )

    # Make sure the output (and plan) directories exist so cp/merge don't fail.
    output_dir.mkdir(parents=True, exist_ok=True)
    if plan.parent != Path(""):
        plan.parent.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(
        [
            {c: e[c] if c != "inputs" else ",".join(e["inputs"]) for c in PLAN_COLUMNS}
            for e in entries
        ],
        columns=PLAN_COLUMNS,
    ).to_csv(plan, sep="\t", index=False)

    n_merge = sum(1 for e in entries if e["action"] == "merge")
    n_copy = sum(1 for e in entries if e["action"] == "copy")
    print(
        f"{len(entries)} groups ({n_merge} merge, {n_copy} copy) written to {plan}.",
        file=sys.stderr,
    )

    slurm_params = dict(
        threads=threads,
        account=account,
        partition=partition,
        time=slurm_time,
        mem_gb=mem_gb,
        mail=mail,
    )

    ready = []
    for e in entries:
        if e["status"] != "ready":
            print(
                f"# SKIP {e['group']} [{e['status']}] missing: {e['missing']}",
                file=sys.stderr,
            )
            continue
        ready.append(
            (
                e["group"],
                group_commands(
                    e,
                    threads,
                    mosdepth=mosdepth,
                    mcg=mcg,
                    bamtocpg=bamtocpg,
                    skip_merge=skip_merge,
                ),
            )
        )

    if slurm and unified:
        # All groups' commands in one sequential SLURM job.
        slurm_dir.mkdir(parents=True, exist_ok=True)
        body = "\n\n".join(
            "# " + group + "\n" + "\n".join(" ".join(c) for c in cmds)
            for group, cmds in ready
        )
        script_path = slurm_dir / "merge_all.slurm.sh"
        script_path.write_text(render_slurm_script("merge_all", body, **slurm_params))
        script_path.chmod(0o755)
        print(
            f"# {len(ready)} groups in one unified script {script_path}. ",
            file=sys.stderr,
        )
        print(f"sbatch {script_path}")
        return 0

    if slurm:
        # One standalone SLURM batch script per group (parallel jobs).
        slurm_dir.mkdir(parents=True, exist_ok=True)
        submit_lines = []
        for group, cmds in ready:
            body = "\n".join(" ".join(c) for c in cmds)
            script_path = slurm_dir / f"merge_{group}.slurm.sh"
            script_path.write_text(
                render_slurm_script(f"merge_{group}", body, **slurm_params)
            )
            script_path.chmod(0o755)
            submit_lines.append(f"sbatch {script_path}")
        print(
            f"# {len(submit_lines)} SLURM scripts written to {slurm_dir}/. ",
            file=sys.stderr,
        )
        for line in submit_lines:
            print(line)
        return 0

    # Plain mode: print (and optionally run) the raw commands.
    for _group, cmds in ready:
        for cmd in cmds:
            print(" ".join(cmd))
            if run and cmd[0] != "#":
                if cmd[0] == "samtools" and shutil.which("samtools") is None:
                    raise ValueError("samtools not found on PATH; cannot --run")
                subprocess.run(cmd, check=True)

    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
