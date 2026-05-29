#!/usr/bin/env bash
# Test for workflow/scripts/make-manifest.py.
# Requires a python with defopt + pandas (the workflow conda env, env.yml).
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
script="$repo/workflow/scripts/make-manifest.py"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# 1) Golden test: generated outputs match the committed expected files.
python "$script" "$here/input.tsv" \
    --manifest "$tmp/manifest.tbl" \
    --decode "$tmp/decode.tsv"
diff "$here/expected_manifest.tbl" "$tmp/manifest.tbl"
diff "$here/expected_decode.tsv" "$tmp/decode.tsv"
echo "PASS: manifest + decode match golden files"

# 2) Unsupported platform (Ultima 'M') must raise.
printf '/x/SMHT005-3A-XX-M59-M001-uwsc-SMAFIZZZ-ug100_GRCh38.aligned.sorted.cram\t/r\t/d\t/o\n' \
    > "$tmp/bad.tsv"
if python "$script" "$tmp/bad.tsv" --manifest "$tmp/x.tbl" --decode "$tmp/x.tsv" 2>/dev/null; then
    echo "FAIL: unsupported platform did not raise" >&2
    exit 1
fi
echo "PASS: unsupported platform raises"
