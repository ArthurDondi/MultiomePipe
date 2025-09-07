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
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to input .tsv.gz fragments file")
    parser.add_argument("--output", required=True, help="Path to output subsampled .tsv.gz file")
    parser.add_argument("--fraction", type=float, default=0.1, help="Fraction of lines to keep")
    parser.add_argument("--cell_data", required=True, help="File with barcodes to keep (one per line)")
    parser.add_argument("--sample", type = str, required=True)
    args = parser.parse_args()

    cell_data = pd.read_csv(args.cell_data, index_col=0)
    cell_data = cell_data[cell_data['sample']==args.sample]
    keep_barcodes = set(cell_data['barcodes'])

    # Step 1: Subsample while preserving top comment lines
    with pysam.BGZFile(args.input, "r") as f_in, pysam.BGZFile(args.output, "w") as f_out:
        # Preserve top comment lines
        for bline in f_in:
            line = bline.decode("utf-8").rstrip("\r\n")
            if line.startswith("#"):
                f_out.write((line + "\n").encode("utf-8"))
                continue
            else:
                parts = line.split("\t")
                barcode = parts[3] 
                if barcode in keep_barcodes:
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
