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


wildcard_constraints:
    tissue="|".join(re.escape(t) for t in GROUPS) or r"(?!.*)",


# ---- targets ----------------------------------------------------------------
TARGETS = []
for _t in GROUPS:
    TARGETS.append(merged_cram(_t))
    TARGETS.append(mosdepth_summary(_t))
    if IS_PACBIO:
        TARGETS.append(mcg_bed(_t))
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


rule cleanup_originals:
    """Opt-in: after the merged CRAM is realized + indexed, delete the per-sample
    results/ CRAMs that fed it. Sentinel-marked so re-runs are no-ops."""
    input:
        cram=lambda wc: merged_cram(wc.tissue),
        crai=lambda wc: merged_cram(wc.tissue) + ".crai",
    output:
        sentinel=touch(os.path.join(OUTPUT_DIR, ".{tissue}.cleaned")),
    params:
        originals=lambda wc: " ".join(
            os.path.join(RESULTS_DIR, f"{s}-{TOKEN}_DSA.aligned.sorted.cram")
            for s in GROUPS[wc.tissue]["samples"]
        ),
    shell:
        r"""
        test -s {input.cram} && test -s {input.crai}
        for f in {params.originals}; do rm -f "$f" "$f".crai; done
        """
