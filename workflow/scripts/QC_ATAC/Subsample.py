#!/usr/bin/env python3
"""
Subsample a fragments file, preserve comments, and create a tabix index.

Usage:
python SubsampleFragments.py \
    --input ../../../input_pbmc_unsorted_3k/ATAC/pbmc_unsorted_3k.atac_fragments.tsv.gz \
    --output pbmc_unsorted_3k.atac_fragments.subsampled.tsv.gz \
    --fraction 0.1
"""

import argparse
import random
import pysam

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to input .tsv.gz fragments file")
    parser.add_argument("--output", required=True, help="Path to output subsampled .tsv.gz file")
    parser.add_argument("--fraction", type=float, default=0.1, help="Fraction of lines to keep")
    args = parser.parse_args()

    # Step 1: Subsample while preserving top comment lines
    with pysam.BGZFile(args.input, "r") as f_in, pysam.BGZFile(args.output, "w") as f_out:
        # Preserve top comment lines
        for bline in f_in:
            line = bline.decode("utf-8").rstrip("\r\n")
            if line.startswith("#"):
                f_out.write((line + "\n").encode("utf-8"))
            else:
                if random.random() < args.fraction:
                    f_out.write((line + "\n").encode("utf-8"))
                break

        # Subsample remaining lines
        for bline in f_in:
            line = bline.decode("utf-8").rstrip("\r\n")
            if random.random() < args.fraction:
                f_out.write((line + "\n").encode("utf-8"))

    # Step 2: Create tabix index
    pysam.tabix_index(
        args.output,
        seq_col=0,        # chromosome column
        start_col=1,      # start column
        end_col=2,        # end column
        meta_char="#",    # comment lines
        force=True,
        zerobased=True
    )

if __name__ == "__main__":
    main()
