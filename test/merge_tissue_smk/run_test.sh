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

echo "ALL MERGE-TISSUE TESTS PASSED"
