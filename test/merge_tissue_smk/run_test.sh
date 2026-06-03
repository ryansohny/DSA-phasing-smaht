#!/usr/bin/env bash
# Dry-run tests for workflow/merge-tissue.smk. Requires `pixi run snakemake`.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
smk="$repo/workflow/merge-tissue.smk"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

R="$tmp/results"
O="$tmp/out"
mkdir -p "$R" "$O"

H="sample\tdonor\tdemographic\ttissue\tpreservation\ttissue_alias\tanalyte\tplatform\tassay\tcenter\tfile_id\taligner\treference\tdsa\tcram"

# Two 3A blood cells (-> merge) + one 3I liver cell (-> copy).
make_decode() {  # $1 = platform string, $2 = outfile
    {
        printf "$H\n"
        printf 'SMHT005-3A-XX-M59-B001-broad-A1\tSMHT005\tM59\tBlood\tSnap Frozen\tblood\tXX\t%s\tWGS\tbroad\tA1\tpbmm2\t/r.fa\t/dsa.fa\t/i1.cram\n' "$1"
        printf 'SMHT005-3A-XX-M59-B001-broad-A2\tSMHT005\tM59\tBlood\tSnap Frozen\tblood\tXX\t%s\tWGS\tbroad\tA2\tpbmm2\t/r.fa\t/dsa.fa\t/i2.cram\n' "$1"
        printf 'SMHT005-3I-001B3-M59-B004-uwsc-I1\tSMHT005\tM59\tLiver\tSnap Frozen\tliver\t001B3\t%s\tFiber-seq\tuwsc\tI1\tpbmm2\t/r.fa\t/dsa.fa\t/i3.cram\n' "$1"
    } >"$2"
}

stub() {  # create per-sample result CRAM stubs
    for s in SMHT005-3A-XX-M59-B001-broad-A1 SMHT005-3A-XX-M59-B001-broad-A2 \
             SMHT005-3I-001B3-M59-B004-uwsc-I1; do
        touch "$R/$s-minimap2_2.31_DSA.aligned.sorted.cram" \
              "$R/$s-minimap2_2.31_DSA.aligned.sorted.cram.crai"
    done
}

dry() {  # run a dry-run; echoes combined output
    pixi run snakemake -n --snakefile "$smk" \
        --config decode="$1" output_dir="$O" results_dir="$R" ${2:-} 2>&1
}

# --- Case 1: PacBio plans merge + mosdepth + mCG -----------------------------
make_decode "PacBio_Revio" "$tmp/pb.tsv"; stub
out=$(dry "$tmp/pb.tsv")
grep -q 'merge_tissue'    <<<"$out" || { echo "FAIL: pacbio no merge_tissue"  >&2; echo "$out" >&2; exit 1; }
grep -q 'mosdepth_tissue' <<<"$out" || { echo "FAIL: pacbio no mosdepth"      >&2; echo "$out" >&2; exit 1; }
grep -q 'mcg_tissue'      <<<"$out" || { echo "FAIL: pacbio no mcg_tissue"    >&2; echo "$out" >&2; exit 1; }
grep -q 'SMHT005-blood-minimap2_2.31_DSA.aligned.sorted.cram' <<<"$out" || { echo "FAIL: alias name" >&2; exit 1; }
echo "PASS: PacBio plans merge + mosdepth + mCG"

# --- Case 1b: Fiber-seq tissue gets a separate Fiber-seq/ merge; WGS does not -
grep -q 'merge_tissue_fiberseq' <<<"$out" || { echo "FAIL: no merge_tissue_fiberseq for fiber tissue" >&2; echo "$out" >&2; exit 1; }
grep -q 'Fiber-seq/SMHT005-liver-minimap2_2.31_DSA.aligned.sorted.cram' <<<"$out" || { echo "FAIL: liver (fiber) missing Fiber-seq/ output" >&2; echo "$out" >&2; exit 1; }
grep -q 'Fiber-seq/SMHT005-blood' <<<"$out" && { echo "FAIL: blood (WGS) must NOT get a Fiber-seq/ output" >&2; echo "$out" >&2; exit 1; }
echo "PASS: Fiber-seq tissue gets a separate Fiber-seq/ merge; WGS tissue does not"

# --- Case 2: ONT plans merge + mosdepth but NOT mCG --------------------------
rm -f "$O"/*.cram "$O"/*.crai 2>/dev/null || true
make_decode "ONT_PromethION24" "$tmp/ont.tsv"
out=$(dry "$tmp/ont.tsv")
grep -q 'mosdepth_tissue' <<<"$out" || { echo "FAIL: ont no mosdepth" >&2; echo "$out" >&2; exit 1; }
grep -q 'mcg_tissue'      <<<"$out" && { echo "FAIL: ont must NOT plan mcg_tissue" >&2; echo "$out" >&2; exit 1; }
echo "PASS: ONT plans merge + mosdepth, no mCG"

# --- Case 3: idempotency — merged present + inputs deleted => nothing to do ---
make_decode "PacBio_Revio" "$tmp/pb.tsv"
mkdir -p "$O/Depth" "$O/mCG"
for t in blood liver; do
    touch "$O/SMHT005-$t-minimap2_2.31_DSA.aligned.sorted.cram" \
          "$O/SMHT005-$t-minimap2_2.31_DSA.aligned.sorted.cram.crai" \
          "$O/Depth/SMHT005-$t-minimap2_2.31_DSA.aligned.sorted.mosdepth.summary.txt" \
          "$O/mCG/SMHT005-$t-minimap2_2.31_DSA.aligned.sorted.combined.bed.gz"
done
# liver (3I) is Fiber-seq, so its Fiber-seq/ merge must also exist for the no-op
mkdir -p "$O/Fiber-seq"
touch "$O/Fiber-seq/SMHT005-liver-minimap2_2.31_DSA.aligned.sorted.cram" \
      "$O/Fiber-seq/SMHT005-liver-minimap2_2.31_DSA.aligned.sorted.cram.crai"
rm -f "$R"/*.cram "$R"/*.crai
out=$(dry "$tmp/pb.tsv")
grep -qiE 'nothing to be done|nothing to do' <<<"$out" || { echo "FAIL: re-run not a no-op" >&2; echo "$out" >&2; exit 1; }
echo "PASS: idempotent re-run after deletion (no MissingInputException)"

# --- Case 4: delete_originals wires cleanup_originals -------------------------
stub  # restore inputs
rm -f "$O"/*.cram "$O"/*.crai "$O"/Depth/* "$O"/mCG/* 2>/dev/null || true
out=$(dry "$tmp/pb.tsv" "delete_originals=true")
grep -q 'cleanup_originals' <<<"$out" || { echo "FAIL: delete_originals did not plan cleanup" >&2; echo "$out" >&2; exit 1; }
echo "PASS: delete_originals plans cleanup_originals"

# --- Case 5: Fiber-seq split details (multi-fiber merge, single-fiber copy, overlap)
Rf="$tmp/rfs"; Of="$tmp/ofs"; mkdir -p "$Rf" "$Of"
{
    printf "$H\n"
    printf 'SMHT005-3A-XX-M59-B001-broad-A1\tSMHT005\tM59\tBlood\tSnap Frozen\tblood\tXX\tPacBio_Revio\tWGS\tbroad\tA1\tpbmm2\t/r.fa\t/dsa.fa\t/i.cram\n'
    printf 'SMHT005-3A-XX-M59-B001-broad-A2\tSMHT005\tM59\tBlood\tSnap Frozen\tblood\tXX\tPacBio_Revio\tWGS\tbroad\tA2\tpbmm2\t/r.fa\t/dsa.fa\t/i.cram\n'
    printf 'SMHT005-3I-001B3-M59-B004-uwsc-I1\tSMHT005\tM59\tLiver\tSnap Frozen\tliver\t001B3\tPacBio_Revio\tFiber-seq\tuwsc\tI1\tpbmm2\t/r.fa\t/dsa.fa\t/i.cram\n'
    printf 'SMHT005-3I-001B3-M59-B004-uwsc-I2\tSMHT005\tM59\tLiver\tSnap Frozen\tliver\t001B3\tPacBio_Revio\tFiber-seq\tuwsc\tI2\tpbmm2\t/r.fa\t/dsa.fa\t/i.cram\n'
    printf 'SMHT005-3C-001A1-M59-B001-broad-C1\tSMHT005\tM59\tKidney\tSnap Frozen\tkidney\t001A1\tPacBio_Revio\tWGS\tbroad\tC1\tpbmm2\t/r.fa\t/dsa.fa\t/i.cram\n'
    printf 'SMHT005-3C-001A2-M59-B004-uwsc-C2\tSMHT005\tM59\tKidney\tSnap Frozen\tkidney\t001A2\tPacBio_Revio\tFiber-seq\tuwsc\tC2\tpbmm2\t/r.fa\t/dsa.fa\t/i.cram\n'
} >"$tmp/fs.tsv"
awk -F'\t' 'NR>1{print $1}' "$tmp/fs.tsv" | while read -r s; do
    touch "$Rf/$s-minimap2_2.31_DSA.aligned.sorted.cram" \
          "$Rf/$s-minimap2_2.31_DSA.aligned.sorted.cram.crai"
done
out=$(pixi run snakemake -n -p --snakefile "$smk" \
    --config decode="$tmp/fs.tsv" output_dir="$Of" results_dir="$Rf" 2>&1)
# liver: 2 fiber cells -> Fiber-seq merge is a samtools merge
grep -E 'samtools merge .*Fiber-seq/SMHT005-liver-minimap2_2.31' <<<"$out" >/dev/null || { echo "FAIL: liver fiber (2 cells) should be a samtools merge into Fiber-seq/" >&2; echo "$out" >&2; exit 1; }
# kidney: 1 fiber cell -> Fiber-seq merge is a cp
grep -E 'cp .*Fiber-seq/SMHT005-kidney-minimap2_2.31' <<<"$out" >/dev/null || { echo "FAIL: kidney fiber (1 cell) should be a cp into Fiber-seq/" >&2; echo "$out" >&2; exit 1; }
# blood (WGS only): no Fiber-seq output
grep -q 'Fiber-seq/SMHT005-blood' <<<"$out" && { echo "FAIL: blood (WGS) must have no Fiber-seq/ output" >&2; echo "$out" >&2; exit 1; }
# overlap: the main kidney merge (output_dir, not Fiber-seq) includes the fiber cell C2
grep -E 'samtools merge .*-o [^ ]*/SMHT005-kidney-minimap2_2.31[^ ]*\.cram .*uwsc-C2' <<<"$out" >/dev/null || { echo "FAIL: main kidney merge should include fiber cell C2 (overlap)" >&2; echo "$out" >&2; exit 1; }
# mosdepth / mCG command lines never operate on a Fiber-seq path
# (anchor on the command flags so Snakemake's "reason:" metadata lines, which list
# both *.mosdepth.summary.txt and Fiber-seq/*.cram, don't false-positive).
grep -E 'mosdepth --threads .*Fiber-seq|aligned_bam_to_cpg_scores --bam .*Fiber-seq' <<<"$out" && { echo "FAIL: mosdepth/mCG must not touch Fiber-seq outputs" >&2; echo "$out" >&2; exit 1; }
echo "PASS: Fiber-seq split (multi-fiber merge, single-fiber copy, overlap, no mosdepth/mCG on fiber)"

echo "ALL MERGE-TISSUE TESTS PASSED"
