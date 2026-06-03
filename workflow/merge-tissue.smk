# Stage-2 workflow: per-tissue merge + mosdepth (+ PacBio mCG) of DSA-phasing
# per-sample CRAMs. Standalone — invoke with its own --configfile, NOT include:d.
import os
import re
import sys
import importlib.util

from snakemake.utils import min_version

min_version("9.0.0")

# ---- config -----------------------------------------------------------------
DECODE = config.get("decode")
if not DECODE:
    raise ValueError("config 'decode' is required (path to the make-manifest decode.tsv)")
OUTPUT_DIR = config.get("output_dir")
if not OUTPUT_DIR:
    raise ValueError("config 'output_dir' is required (destination for merged CRAMs)")
OUTPUT_DIR = os.path.abspath(OUTPUT_DIR)
RESULTS_DIR = config.get("results_dir", "results")
DELETE_ORIGINALS = bool(config.get("delete_originals", False))
THREADS = int(config.get("threads", 39))
MOSDEPTH_THREADS = min(4, THREADS)

DEFAULT_ENV = "envs/env.yml"

# ---- import the standalone tool (hyphenated filename) -----------------------
_script = os.path.join(workflow.basedir, "scripts", "merge-by-tissue-pacbio-n-extract-mcg.py")
_spec = importlib.util.spec_from_file_location("merge_by_tissue", _script)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
build_plan = _mod.build_plan
merge_command = _mod.merge_command
mosdepth_commands = _mod.mosdepth_commands
mcg_commands = _mod.mcg_commands
BAMTOCPG = config.get("bamtocpg", _mod.DEFAULT_BAMTOCPG)


# ---- deterministic aligner token (does NOT depend on inputs existing) -------
def _minimap2_version():
    path = os.path.join(workflow.basedir, "envs", "env.yml")
    with open(path) as fh:
        m = re.search(r"minimap2==([0-9][0-9A-Za-z.\-]*)", fh.read())
    if not m:
        raise ValueError(f"No pinned minimap2 version in {path}")
    return m.group(1)


TOKEN = f"minimap2_{_minimap2_version()}"

# ---- grouping ---------------------------------------------------------------
ENTRIES = build_plan(DECODE, RESULTS_DIR, OUTPUT_DIR, ["tissue_code"])
GROUPS = {e["name"]: e for e in ENTRIES}

import pandas as pd

_df = pd.read_csv(DECODE, sep="\t", dtype=str)
_donors = set(_df["donor"])
if len(_donors) != 1:
    raise ValueError(f"decode spans multiple donors {sorted(_donors)}; one donor per run")
DONOR = _donors.pop()
IS_PACBIO = all(p.startswith("PacBio") for p in _df["platform"])

# Fiber-seq cells (decode 'assay' == "Fiber-seq"), grouped per tissue. Fiber cells
# are ALSO merged separately into output_dir/Fiber-seq/ (kept, no mosdepth/mCG),
# in addition to being part of the main all-cells tissue merge.
_assay = dict(zip(_df["sample"], _df["assay"]))


def _is_fiber(sample):
    return _assay.get(sample, "") == "Fiber-seq"


FIBER_GROUPS = {t: [s for s in GROUPS[t]["samples"] if _is_fiber(s)] for t in GROUPS}
FIBER_GROUPS = {t: ss for t, ss in FIBER_GROUPS.items() if ss}


def merged_cram(tissue):
    return os.path.join(OUTPUT_DIR, f"{DONOR}-{tissue}-{TOKEN}_DSA.aligned.sorted.cram")


def mosdepth_summary(tissue):
    base = os.path.basename(merged_cram(tissue))[: -len(".cram")]
    return os.path.join(OUTPUT_DIR, "Depth", base + ".mosdepth.summary.txt")


def mcg_bed(tissue):
    base = os.path.basename(merged_cram(tissue))[: -len(".cram")]
    return os.path.join(OUTPUT_DIR, "mCG", base + ".combined.bed.gz")


def cleaned_sentinel(tissue):
    return os.path.join(OUTPUT_DIR, f".{tissue}.cleaned")


def fiber_cram(tissue):
    return os.path.join(
        OUTPUT_DIR, "Fiber-seq", f"{DONOR}-{tissue}-{TOKEN}_DSA.aligned.sorted.cram"
    )


wildcard_constraints:
    tissue="|".join(re.escape(t) for t in GROUPS) or r"(?!.*)",


# ---- targets ----------------------------------------------------------------
TARGETS = []
for _t in GROUPS:
    TARGETS.append(merged_cram(_t))
    TARGETS.append(mosdepth_summary(_t))
    if IS_PACBIO:
        TARGETS.append(mcg_bed(_t))
    if _t in FIBER_GROUPS:
        TARGETS.append(fiber_cram(_t))
    if DELETE_ORIGINALS:
        TARGETS.append(cleaned_sentinel(_t))

print(
    f"[merge-tissue] donor={DONOR} platform={'PacBio' if IS_PACBIO else 'ONT'} "
    f"groups={len(GROUPS)} token={TOKEN} delete_originals={DELETE_ORIGINALS}",
    file=sys.stderr,
)


rule all:
    input:
        TARGETS


def merge_inputs(wc):
    """Per-sample CRAMs to merge — but empty once the merged CRAM already exists,
    so a re-run after the inputs were deleted is a clean no-op instead of a
    MissingInputException."""
    if os.path.exists(merged_cram(wc.tissue)):
        return []
    return [
        os.path.join(RESULTS_DIR, f"{s}-{TOKEN}_DSA.aligned.sorted.cram")
        for s in GROUPS[wc.tissue]["samples"]
    ]


def _merge_cmd(wc, input, output):
    crams = list(input)
    if not crams:
        return "true"  # merged CRAM already present; rule will not execute
    if GROUPS[wc.tissue]["action"] == "copy":
        src = crams[0]
        return f"cp {src} {output.cram} && cp {src}.crai {output.crai}"
    cmd = merge_command({"output": output.cram, "inputs": crams}, THREADS)
    return " ".join(cmd)


rule merge_tissue:
    input:
        crams=merge_inputs,
    output:
        cram=os.path.join(OUTPUT_DIR, f"{DONOR}-{{tissue}}-{TOKEN}_DSA.aligned.sorted.cram"),
        crai=os.path.join(OUTPUT_DIR, f"{DONOR}-{{tissue}}-{TOKEN}_DSA.aligned.sorted.cram.crai"),
    params:
        cmd=lambda wc, input, output: _merge_cmd(wc, input, output),
    conda:
        DEFAULT_ENV
    threads: THREADS
    resources:
        runtime=24 * 60,
        mem_mb=THREADS * 1024,
    shell:
        "{params.cmd}"


rule mosdepth_tissue:
    input:
        cram=lambda wc: merged_cram(wc.tissue),
    output:
        summary=os.path.join(
            OUTPUT_DIR, "Depth", f"{DONOR}-{{tissue}}-{TOKEN}_DSA.aligned.sorted.mosdepth.summary.txt"
        ),
    params:
        cmd=lambda wc: " && ".join(
            " ".join(c)
            for c in mosdepth_commands(merged_cram(wc.tissue), GROUPS[wc.tissue]["dsa"], MOSDEPTH_THREADS)
        ),
    conda:
        DEFAULT_ENV
    threads: MOSDEPTH_THREADS
    resources:
        runtime=12 * 60,
        mem_mb=8 * 1024,
    shell:
        "{params.cmd}"


rule mcg_tissue:
    input:
        cram=lambda wc: merged_cram(wc.tissue),
    output:
        bed=os.path.join(
            OUTPUT_DIR, "mCG", f"{DONOR}-{{tissue}}-{TOKEN}_DSA.aligned.sorted.combined.bed.gz"
        ),
    params:
        cmd=lambda wc: " && ".join(
            " ".join(c)
            for c in mcg_commands(merged_cram(wc.tissue), GROUPS[wc.tissue]["dsa"], THREADS, BAMTOCPG)
        ),
    conda:
        DEFAULT_ENV
    threads: THREADS
    resources:
        runtime=24 * 60,
        mem_mb=THREADS * 1024,
    shell:
        "{params.cmd}"


def fiber_merge_inputs(wc):
    """Fiber-seq per-sample CRAMs for the tissue; empty once the Fiber-seq merge
    already exists (idempotent re-run after deletion)."""
    if os.path.exists(fiber_cram(wc.tissue)):
        return []
    return [
        os.path.join(RESULTS_DIR, f"{s}-{TOKEN}_DSA.aligned.sorted.cram")
        for s in FIBER_GROUPS[wc.tissue]
    ]


def _fiber_merge_cmd(wc, input, output):
    crams = list(input)
    if not crams:
        return "true"  # fiber merge already present; rule will not execute
    if len(crams) == 1:
        src = crams[0]
        return f"cp {src} {output.cram} && cp {src}.crai {output.crai}"
    cmd = merge_command({"output": output.cram, "inputs": crams}, THREADS)
    return " ".join(cmd)


rule merge_tissue_fiberseq:
    """Merge the tissue's Fiber-seq cells into output_dir/Fiber-seq/ (kept; no
    mosdepth/mCG). Fiber cells also appear in the main all-cells tissue merge."""
    input:
        crams=fiber_merge_inputs,
    output:
        cram=os.path.join(
            OUTPUT_DIR, "Fiber-seq", f"{DONOR}-{{tissue}}-{TOKEN}_DSA.aligned.sorted.cram"
        ),
        crai=os.path.join(
            OUTPUT_DIR,
            "Fiber-seq",
            f"{DONOR}-{{tissue}}-{TOKEN}_DSA.aligned.sorted.cram.crai",
        ),
    wildcard_constraints:
        tissue="|".join(re.escape(t) for t in FIBER_GROUPS) or r"(?!.*)",
    params:
        cmd=lambda wc, input, output: _fiber_merge_cmd(wc, input, output),
    conda:
        DEFAULT_ENV
    threads: THREADS
    resources:
        runtime=24 * 60,
        mem_mb=THREADS * 1024,
    shell:
        "{params.cmd}"


def cleanup_deps(wc):
    """Merged outputs that must exist (and be verified) before a tissue's
    per-sample results/ CRAMs are deleted: the main all-cells merge, plus the
    Fiber-seq merge when the tissue has fiber cells (fiber cells feed both)."""
    deps = [merged_cram(wc.tissue), merged_cram(wc.tissue) + ".crai"]
    if wc.tissue in FIBER_GROUPS:
        deps += [fiber_cram(wc.tissue), fiber_cram(wc.tissue) + ".crai"]
    return deps


rule cleanup_originals:
    """Opt-in: after every merge that consumes a tissue's cells is realized +
    indexed, delete ALL of that tissue's per-sample results/ CRAMs (WGS and
    Fiber-seq). Sentinel-marked so re-runs are no-ops."""
    input:
        merged=cleanup_deps,
    output:
        sentinel=touch(os.path.join(OUTPUT_DIR, ".{tissue}.cleaned")),
    params:
        originals=lambda wc: " ".join(
            os.path.join(RESULTS_DIR, f"{s}-{TOKEN}_DSA.aligned.sorted.cram")
            for s in GROUPS[wc.tissue]["samples"]
        ),
    shell:
        r"""
        for f in {input.merged}; do test -s "$f"; done
        for f in {params.originals}; do rm -f "$f" "$f".crai; done
        """
