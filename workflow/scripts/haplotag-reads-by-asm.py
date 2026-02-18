#!/usr/bin/env python
import logging
import sys
from pathlib import Path
from typing import Optional

import defopt
import pysam
from tqdm import tqdm
import pandas as pd


def run(
    bam, obam, hap1_tag, hap2_tag, min_mapq, reset_mapq=None, reset_mapq_before=False
):
    read_assignments = []
    for rec in tqdm(bam.fetch(until_eof=True)):
        # save original MAPQ in om tag if we will reset it
        if reset_mapq is not None:
            rec.set_tag("om", rec.mapping_quality)

        # optionally reset MAPQ before haplotag assignment
        if reset_mapq is not None and reset_mapq_before and not rec.is_unmapped:
            rec.mapping_quality = reset_mapq

        tag_value = None
        if rec.is_unmapped:
            tag_value = None
        else:
            matches_h1 = hap1_tag in rec.reference_name
            matches_h2 = hap2_tag in rec.reference_name
            matches_both = matches_h1 and matches_h2
            if matches_both:
                tag_value = None
            elif matches_h1:
                tag_value = 1
            elif matches_h2:
                tag_value = 2
            else:
                tag_value = None
        # set the oh tag first
        rec.set_tag("oh", tag_value)
        # if we have a low mapq reset just the HP tag
        if rec.mapping_quality < min_mapq:
            tag_value = None
        rec.set_tag("HP", tag_value)

        # optionally reset MAPQ after haplotag assignment
        if reset_mapq is not None and not reset_mapq_before and not rec.is_unmapped:
            rec.mapping_quality = reset_mapq

        is_primary = not rec.is_secondary and not rec.is_supplementary
        if is_primary and tag_value is not None:
            read_assignments.append((tag_value, rec.query_name))

        obam.write(rec)
    return read_assignments


def main(
    infile: Optional[Path] = None,
    outfile: Optional[Path] = None,
    outfile2: Optional[Path] = None,
    *,
    hap1_tag: str = "haplotype1",
    hap2_tag: str = "haplotype2",
    min_mapq: int = 1,
    reset_mapq: Optional[int] = None,
    reset_mapq_before: bool = False,
    threads: int = 8,
    verbose: int = 0,
):
    """
    Assign haplotype tags to reads based on DSA contig names.

    Reads aligned to a DSA are tagged with HP (haplotype) based on whether
    the reference contig name contains the hap1 or hap2 tag string. An
    additional "oh" (original haplotype) tag preserves the assignment before
    MAPQ filtering, so low-MAPQ reads lose HP but retain "oh" for debugging.

    When reset_mapq is set, original MAPQ is saved in the "om" tag and MAPQ
    is overwritten to the specified value. Timing is controlled by
    reset_mapq_before: before haplotagging (min_mapq filter uses the new
    value) or after (min_mapq filter uses the original value).

    :param infile: Input BAM file, stdin by default
    :param outfile: Output BAM file, stdout by default
    :param outfile2: Output TSV file for read-to-haplotype assignments
    :param hap1_tag: Substring in contig names identifying haplotype 1
    :param hap2_tag: Substring in contig names identifying haplotype 2
    :param min_mapq: Minimum mapping quality; reads below this keep "oh" but lose HP
    :param reset_mapq: Reset MAPQ to this value; original saved in "om" tag
    :param reset_mapq_before: If true, reset MAPQ before haplotag assignment
    :param threads: Number of threads for BAM I/O
    :param verbose: Set the logging level of the function
    """
    logger = logging.getLogger()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logger.setLevel(log_level)

    if infile is None or infile == "-":
        infile = sys.stdin
    if outfile is None or outfile == "-":
        outfile = sys.stdout

    bam = pysam.AlignmentFile(infile, "rb", threads=threads)
    obam = pysam.AlignmentFile(outfile, "wb", template=bam, threads=threads)
    read_assignments = run(
        bam,
        obam,
        hap1_tag,
        hap2_tag,
        min_mapq,
        reset_mapq=reset_mapq,
        reset_mapq_before=reset_mapq_before,
    )

    pd.DataFrame(read_assignments, columns=["haplotype", "read_name"]).to_csv(
        outfile2, index=False, sep="\t"
    )

    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
