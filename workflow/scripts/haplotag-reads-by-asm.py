#!/usr/bin/env python
import logging
import sys
from pathlib import Path
from typing import Optional

import defopt
import pysam
from tqdm import tqdm
import pandas as pd


def run(bam, obam, hap1_tag, hap2_tag, min_mapq):
    read_assignments = []
    for rec in tqdm(bam.fetch(until_eof=True)):
        is_primary = not rec.is_secondary and not rec.is_supplementary
        matches_h1 = hap1_tag in rec.reference_name
        matches_h2 = hap2_tag in rec.reference_name
        matches_both = matches_h1 and matches_h2
        
        tag_value = None
        if rec.is_unmapped or matches_both: 
            tag_value = None
        else:
            if matches_h1:
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
    threads: int = 8,
    verbose: int = 0,
):
    """
    Author Mitchell R. Vollger

    :param infile: Input file, stdin by default
    :param outfile: Output file, stdout by default
    :param outfile2: Output fil
    :param hap1_tag: Haplotype 1 contig tag
    :param hap2_tag: Haplotype 2 contig tag
    :param min_mapq: Minimum mapping quality to consider
    :param verbose: Set the logging level of the function
    """
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
    )

    pd.DataFrame(read_assignments, columns=["haplotype", "read_name"]).to_csv(
        outfile2, index=False, sep="\t"
    )

    logger = logging.getLogger()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logger.setLevel(log_level)

    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
