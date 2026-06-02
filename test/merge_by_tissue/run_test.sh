#!/usr/bin/env bash
# Test for workflow/scripts/merge-by-tissue-pacbio-n-extract-mcg.py.
# Requires a python with defopt + pandas (the workflow conda env, env.yml).
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
script="$repo/workflow/scripts/merge-by-tissue-pacbio-n-extract-mcg.py"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

R="$tmp/results"
O="$tmp/merged"
mkdir -p "$R" "$O"

# decode table: two 3A cells (-> merge) + one 3I cell (-> copy)
H="sample	donor	demographic	tissue	preservation	tissue_alias	analyte	platform	assay	center	file_id	aligner	reference	dsa	cram"
{
    printf '%s\n' "$H"
    printf 'SMHT005-3A-XX-M59-B001-broad-SMAFICOSSIZL\tSMHT005\tM59\tBlood\tSnap Frozen\tblood\tXX\tPacBio_Revio\tWGS\tbroad\tSMAFICOSSIZL\tpbmm2_1.13.0\t/ref.fa\t/dsa.fa\t/in1.cram\n'
    printf 'SMHT005-3A-XX-M59-B001-broad-SMAFIO1QBEN4\tSMHT005\tM59\tBlood\tSnap Frozen\tblood\tXX\tPacBio_Revio\tWGS\tbroad\tSMAFIO1QBEN4\tpbmm2_1.13.0\t/ref.fa\t/dsa.fa\t/in2.cram\n'
    printf 'SMHT005-3I-001B3-M59-B004-uwsc-SMAFI8EOWRG5\tSMHT005\tM59\tLiver\tSnap Frozen\tliver\t001B3\tPacBio_Revio\tFiber-seq\tuwsc\tSMAFI8EOWRG5\tpbmm2_1.13.0\t/ref.fa\t/dsa.fa\t/in3.cram\n'
} >"$tmp/decode.tsv"

for s in SMHT005-3A-XX-M59-B001-broad-SMAFICOSSIZL \
    SMHT005-3A-XX-M59-B001-broad-SMAFIO1QBEN4 \
    SMHT005-3I-001B3-M59-B004-uwsc-SMAFI8EOWRG5; do
    touch "$R/$s-minimap2_2.31_DSA.aligned.sorted.cram" \
        "$R/$s-minimap2_2.31_DSA.aligned.sorted.cram.crai"
done

python "$script" "$tmp/decode.tsv" \
    --results-dir "$R" --output-dir "$O" --group-by tissue_code \
    --plan "$tmp/plan.tsv" >/dev/null 2>&1

# Expect: 3A merge (2 inputs, ready), 3I copy (1 input, ready)
grep -qP '^3A\tmerge\tready\t2\t' "$tmp/plan.tsv" || { echo "FAIL: 3A not a 2-input merge" >&2; cat "$tmp/plan.tsv" >&2; exit 1; }
grep -qP '^3I\tcopy\tready\t1\t' "$tmp/plan.tsv" || { echo "FAIL: 3I not a 1-input copy" >&2; cat "$tmp/plan.tsv" >&2; exit 1; }
grep -q 'SMHT005-blood-minimap2_2.31_DSA.aligned.sorted.cram' "$tmp/plan.tsv" || { echo "FAIL: merged name should use tissue alias" >&2; cat "$tmp/plan.tsv" >&2; exit 1; }
echo "PASS: tissue grouping -> merge/copy plan correct"

# SLURM mode: per-group batch scripts following the lab template.
python "$script" "$tmp/decode.tsv" \
    --results-dir "$R" --output-dir "$O" --group-by tissue_code \
    --threads 20 --mem-gb 64 --slurm --slurm-dir "$tmp/slurm" \
    --plan "$tmp/plan2.tsv" >/dev/null 2>&1
s="$tmp/slurm/merge_3A.slurm.sh"
[ -f "$s" ] || { echo "FAIL: SLURM script not written" >&2; exit 1; }
grep -q '^#SBATCH --account=stergachislab' "$s" || { echo "FAIL: missing SBATCH account" >&2; exit 1; }
grep -q '^#SBATCH --cpus-per-task=20' "$s" || { echo "FAIL: cpus-per-task not from --threads" >&2; exit 1; }
grep -q '^#SBATCH --mem=64G' "$s" || { echo "FAIL: mem not from --mem-gb" >&2; exit 1; }
grep -q 'samtools merge' "$s" || { echo "FAIL: merge command missing from script" >&2; exit 1; }
grep -q -- '-O CRAM' "$s" || { echo "FAIL: merge not using -O CRAM" >&2; exit 1; }
grep -q -- '--reference' "$s" && { echo "FAIL: merge should not pass --reference (embed_ref inputs)" >&2; exit 1; }
# mosdepth + mCG steps present by default, on the merged CRAM
grep -q 'mosdepth .*--no-per-base' "$s" || { echo "FAIL: mosdepth step missing" >&2; exit 1; }
grep -q 'aligned_bam_to_cpg_scores .*--hap-tag HP' "$s" || { echo "FAIL: mCG step missing" >&2; exit 1; }
grep -q 'mosdepth --threads 4' "$s" || { echo "FAIL: mosdepth should cap at 4 threads" >&2; exit 1; }
echo "PASS: SLURM scripts follow the template (+ mosdepth + mCG)"

# Toggles: --no-mosdepth --no-mcg leaves only the merge.
python "$script" "$tmp/decode.tsv" \
    --results-dir "$R" --output-dir "$O" --group-by tissue_code \
    --threads 20 --slurm --slurm-dir "$tmp/slurm_noextra" --no-mosdepth --no-mcg \
    --plan "$tmp/plan2b.tsv" >/dev/null 2>&1
sn="$tmp/slurm_noextra/merge_3A.slurm.sh"
grep -q 'samtools merge' "$sn" || { echo "FAIL: merge missing with toggles off" >&2; exit 1; }
grep -q 'mosdepth' "$sn" && { echo "FAIL: --no-mosdepth should drop mosdepth" >&2; exit 1; }
grep -q 'aligned_bam_to_cpg_scores' "$sn" && { echo "FAIL: --no-mcg should drop mCG" >&2; exit 1; }
echo "PASS: --no-mosdepth/--no-mcg toggles work"

# --skip-merge: run only mosdepth/mCG on an already-merged CRAM.
touch "$O/SMHT005-blood-minimap2_2.31_DSA.aligned.sorted.cram" \
    "$O/SMHT005-blood-minimap2_2.31_DSA.aligned.sorted.cram.crai"
python "$script" "$tmp/decode.tsv" \
    --results-dir "$R" --output-dir "$O" --group-by tissue_code \
    --threads 20 --slurm --slurm-dir "$tmp/slurm_skip" --skip-merge \
    --plan "$tmp/plan_skip.tsv" >/dev/null 2>&1
ss="$tmp/slurm_skip/merge_3A.slurm.sh"
grep -qP '^3A\tmerge\tready\t' "$tmp/plan_skip.tsv" || { echo "FAIL: skip-merge 3A not ready (merged CRAM exists)" >&2; cat "$tmp/plan_skip.tsv" >&2; exit 1; }
grep -q '^# samtools merge' "$ss" || { echo "FAIL: --skip-merge should keep merge commented" >&2; exit 1; }
grep -qE '^samtools merge' "$ss" && { echo "FAIL: --skip-merge merge must not be executable" >&2; exit 1; }
grep -q 'mosdepth .*--no-per-base' "$ss" || { echo "FAIL: skip-merge missing mosdepth" >&2; exit 1; }
grep -q 'aligned_bam_to_cpg_scores' "$ss" || { echo "FAIL: skip-merge missing mCG" >&2; exit 1; }
echo "PASS: --skip-merge comments the merge, runs mosdepth/mCG"

# Unified mode: one script with all groups' commands.
python "$script" "$tmp/decode.tsv" \
    --results-dir "$R" --output-dir "$O" --group-by tissue_code \
    --threads 16 --slurm --unified --slurm-dir "$tmp/uni" \
    --plan "$tmp/plan3.tsv" >/dev/null 2>&1
u="$tmp/uni/merge_all.slurm.sh"
[ -f "$u" ] || { echo "FAIL: unified script not written" >&2; exit 1; }
[ "$(ls "$tmp/uni" | wc -l)" -eq 1 ] || { echo "FAIL: unified should write exactly one script" >&2; exit 1; }
grep -q '^#SBATCH --job-name=merge_all' "$u" || { echo "FAIL: unified job-name wrong" >&2; exit 1; }
grep -q 'SMHT005-blood-minimap2_2.31_DSA.aligned.sorted.cram' "$u" || { echo "FAIL: unified missing 3A merge" >&2; exit 1; }
grep -q 'SMHT005-liver-minimap2_2.31_DSA.aligned.sorted.cram' "$u" || { echo "FAIL: unified missing 3I copy" >&2; exit 1; }
echo "PASS: unified SLURM script contains all groups"
