#!/usr/bin/env python
import sys
from pathlib import Path
from typing import Optional

import defopt
import pysam
from tqdm import tqdm
import pandas as pd
import logging

def main(
    infile: Optional[Path] = None,
    assignments: Optional[Path] = None,
    outfile: Optional[Path] = None,
    *,
    threads: int = 8,
    verbose: int = 0,
):
    """
    Author Mitchell R. Vollger

    :param infile: Input file, stdin by default
    :param assignments: haplotype assignments for reads
    :param outfile: Output file, stdout by default
    :param count: Number of times to display the greeting
    :param verbose: Set the logging level of the function
    """
    if infile is None or infile == "-":
        infile = sys.stdin
    if outfile is None or outfile == "-":
        outfile = sys.stdout

    # make the query name the index
    assignments = pd.read_csv(assignments, sep="\t").set_index("read_name")

    bam = pysam.AlignmentFile(infile, "rb", threads=threads)
    o_bam = pysam.AlignmentFile(outfile, "wbu", template=bam, threads=threads)
    for rec in tqdm(bam.fetch(until_eof=True)):
        if rec.query_name in assignments.index:
            hp = assignments.loc[rec.query_name, "haplotype"]
            logging.warn(f"Setting HP tag to {hp} for {rec.query_name}")
            rec.set_tag("HP", hp)
        else:
            rec.set_tag("HP", None)
        o_bam.write(rec)

    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
