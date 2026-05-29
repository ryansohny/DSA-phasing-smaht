#!/usr/bin/env python
# Generated with assist from Claude code (Opus 4.8)
"""Generate a DSA-phasing manifest.tbl (and a decoded metadata table) from a TSV
of SMaHT CRAM rows.

The input TSV has 4 tab-separated columns and no header:

    CRAM    hg38_ref    DSA    output_dir

Each CRAM filename is decoded against the SMaHT lookup tables (assay, tissue,
platform) bundled in the lookups directory. The platform and Fiber-seq status
come from the platform+assay token in the filename (e.g. B004 = PacBio +
Fiber-seq).
"""

import logging
import sys
from pathlib import Path
from typing import Optional

import defopt
import pandas as pd

LOOKUP_DIR = Path(__file__).resolve().parent / "lookups"

# Platform_lookup codes (first char of the platform+assay token) collapsed to
# the three platforms the DSA-phasing pipeline understands.
PLATFORM_CLASS = {
    "A": "illumina",  # Illumina_NovaSeqX
    "B": "pacbio",  # PacBio_Revio
    "C": "illumina",  # Illumina_NovaSeq6000
    "D": "ont",  # ONT_PromethION24
    "E": "ont",  # ONT_PromethION2Solo
    "F": "ont",  # ONT_MiniONMk1B
    "G": "illumina",  # Illumina_HiSeqX
    "K": "illumina",  # Illumina_NextSeq2000
    "L": "pacbio",  # PacBio_SequelIIe
}

# Assay code that marks Fiber-seq (Assay_lookup.tsv).
FIBERSEQ_ASSAY = "004"

MANIFEST_COLUMNS = [
    "sample",
    "dsa",
    "bam",
    "h1_tag",
    "h2_tag",
    "platform",
    "fiber-seq",
]

DECODE_COLUMNS = [
    "sample",
    "donor",
    "demographic",
    "tissue",
    "preservation",
    "tissue_alias",
    "analyte",
    "platform",
    "assay",
    "center",
    "file_id",
    "aligner",
    "reference",
    "dsa",
    "cram",
]


def load_lookups(lookup_dir):
    """Load the assay/tissue/platform lookup tables into dicts keyed by Code."""
    assay = pd.read_csv(lookup_dir / "Assay_lookup.tsv", sep="\t", dtype=str)
    tissue = pd.read_csv(lookup_dir / "Tissue_lookup.tsv", sep="\t", dtype=str)
    platform = pd.read_csv(lookup_dir / "Platform_lookup.tsv", sep="\t", dtype=str)
    return {
        "assay": assay.set_index("Code")["Assay"].to_dict(),
        "tissue": tissue.set_index("Code").to_dict("index"),
        "platform": platform.set_index("Code")["Platform"].to_dict(),
    }


def sample_name_from_cram(cram_path):
    """Strip the trailing `-<aligner>_...cram` segment to get the sample name.

    The aligner token (e.g. `pbmm2_1.13.0_GRCh38.aligned.sorted.cram`) is the
    last dash-delimited segment, so the sample name is everything before it.
    Returns (sample_name, aligner_token).
    """
    base = Path(cram_path).name
    if "-" not in base:
        raise ValueError(f"CRAM name has no dash-delimited fields: {base!r}")
    sample, aligner_segment = base.rsplit("-", 1)
    # aligner is the part before `_GRCh38...`; fall back to the whole segment.
    aligner = aligner_segment.split("_GRCh38")[0]
    return sample, aligner


def decode_sample(sample, lookups):
    """Decode a SMaHT sample name into its constituent metadata fields."""
    fields = sample.split("-")
    if len(fields) != 7:
        raise ValueError(
            f"Expected 7 dash-delimited fields in sample {sample!r}, got {len(fields)}"
        )
    donor, tissue_code, analyte, demographic, plat_assay, center, file_id = fields

    if len(plat_assay) < 4:
        raise ValueError(
            f"Sample {sample!r}: platform+assay token {plat_assay!r} is too short"
        )
    plat_code, assay_code = plat_assay[0], plat_assay[1:]

    if plat_code not in lookups["platform"]:
        raise ValueError(
            f"Sample {sample!r}: unknown platform code {plat_code!r} "
            f"(token {plat_assay!r})"
        )
    if plat_code not in PLATFORM_CLASS:
        raise ValueError(
            f"Sample {sample!r}: platform {lookups['platform'][plat_code]!r} "
            f"(code {plat_code!r}) is not supported by the DSA-phasing pipeline "
            f"(pacbio/ont/illumina only)"
        )
    if assay_code not in lookups["assay"]:
        raise ValueError(
            f"Sample {sample!r}: unknown assay code {assay_code!r} "
            f"(token {plat_assay!r})"
        )
    if tissue_code not in lookups["tissue"]:
        raise ValueError(f"Sample {sample!r}: unknown tissue code {tissue_code!r}")

    tissue_row = lookups["tissue"][tissue_code]
    return {
        "donor": donor,
        "demographic": demographic,
        "tissue": tissue_row["Tissue"],
        "preservation": tissue_row["Preservation"],
        "tissue_alias": tissue_row["Tissue_alias"],
        "analyte": analyte,
        "platform_class": PLATFORM_CLASS[plat_code],
        "platform_full": lookups["platform"][plat_code],
        "assay": lookups["assay"][assay_code],
        "fiber_seq": assay_code == FIBERSEQ_ASSAY,
        "center": center,
        "file_id": file_id,
    }


def build_rows(input_tsv, lookups):
    """Parse the input TSV and decode each row. Returns (manifest, decode) lists."""
    manifest_rows = []
    decode_rows = []
    references = set()

    df = pd.read_csv(input_tsv, sep="\t", header=None, dtype=str)
    if df.shape[1] != 4:
        raise ValueError(
            f"Input {input_tsv} has {df.shape[1]} columns, expected 4 "
            f"(CRAM, hg38_ref, DSA, output_dir)"
        )

    for cram, reference, dsa, _outdir in df.itertuples(index=False, name=None):
        sample, aligner = sample_name_from_cram(cram)
        meta = decode_sample(sample, lookups)
        references.add(reference)

        manifest_rows.append(
            {
                "sample": sample,
                "dsa": dsa,
                "bam": cram,
                "h1_tag": "h1",
                "h2_tag": "h2",
                "platform": meta["platform_class"],
                "fiber-seq": "true" if meta["fiber_seq"] else "false",
            }
        )
        decode_rows.append(
            {
                "sample": sample,
                "donor": meta["donor"],
                "demographic": meta["demographic"],
                "tissue": meta["tissue"],
                "preservation": meta["preservation"],
                "tissue_alias": meta["tissue_alias"],
                "analyte": meta["analyte"],
                "platform": meta["platform_full"],
                "assay": meta["assay"],
                "center": meta["center"],
                "file_id": meta["file_id"],
                "aligner": aligner,
                "reference": reference,
                "dsa": dsa,
                "cram": cram,
            }
        )

    if len(references) > 1:
        raise ValueError(
            f"Input rows reference multiple hg38 references: {sorted(references)}. "
            f"The pipeline `cram_reference` is a single global value."
        )

    return manifest_rows, decode_rows, references.pop() if references else None


def main(
    input_tsv: Path,
    *,
    manifest: Path = Path("manifest.tbl"),
    decode: Path = Path("decode.tsv"),
    lookup_dir: Optional[Path] = None,
    verbose: int = 0,
):
    """
    Generate a DSA-phasing manifest.tbl and a decoded metadata table from a
    4-column SMaHT CRAM TSV (CRAM, hg38_ref, DSA, output_dir).

    :param input_tsv: 4-column tab-separated input (no header)
    :param manifest: Output path for the pipeline manifest.tbl
    :param decode: Output path for the human-readable decoded metadata table
    :param lookup_dir: Directory holding the lookup tables (defaults to bundled)
    :param verbose: Set the logging level of the function
    """
    logging.basicConfig(
        format="[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    )
    logging.getLogger().setLevel(10 * (3 - verbose))

    lookups = load_lookups(lookup_dir or LOOKUP_DIR)
    manifest_rows, decode_rows, reference = build_rows(input_tsv, lookups)

    pd.DataFrame(manifest_rows, columns=MANIFEST_COLUMNS).to_csv(
        manifest, sep="\t", index=False
    )
    pd.DataFrame(decode_rows, columns=DECODE_COLUMNS).to_csv(
        decode, sep="\t", index=False
    )

    print(
        f"Wrote {len(manifest_rows)} samples to {manifest} and {decode}.",
        file=sys.stderr,
    )
    if reference is not None:
        print(
            f"All inputs share one reference. Add to your config:\n"
            f"    cram_reference: {reference}",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
