#!/usr/bin/env python3
r"""
Apply QC filters to a single sample and generate QC plots.

This script:
- Computes barcodes passing QC thresholds (user-defined or automatic)
- Generates sample-level statistics plots
- Generates barcode-level statistics plots (with detailed info)
- Writes the barcodes passing QC to a text file (one barcode per line)

Usage:
python CreateATACCountMatrix.py \
    --sample sampleA \
    --barcodes barcodes_passing_qc.tsv \
    --fragments path/to/fragments.tsv.gz \
    --regions path/to/consensus_peaks.bed \
    --blacklist path/to/blacklist.bed \
    --parquet path/to/sample_metrics.parquet \
    --n_cpu 8 \
"""

import argparse
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import polars as pl
import pickle


def main():
    p = argparse.ArgumentParser(description="Apply QC filters and build cisTopic object.")
    p.add_argument("--sample", type=str, required=True, help="Sample name.")
    p.add_argument("--barcodes", required=True, help="File with barcodes passing QC (one per line).")
    p.add_argument("--fragments", required=True, help="Path to fragments file (.tsv.gz).")
    p.add_argument("--out_cistopic_obj", required=True, help="Path to output pickled cistopic object.")
    p.add_argument("--regions", required=True, help="Path to consensus peaks BED file.")
    p.add_argument("--blacklist", required=True, help="Path to blacklist BED file.")
    p.add_argument("--parquet", required=True, help="Path to metrics parquet file.")
    p.add_argument("--min_fragments_per_region", type=int, default=0, help="Min. fragments in peak to keep peak")
    p.add_argument("--min_cells_per_region", type=int, default=0, help="Min. cells in peak to keep peak")
    p.add_argument("--n_cpu", type=int, default=1, help="Number of CPUs to use.")
    p.add_argument("--split_pattern", type=str, default='-')

    args = p.parse_args()
    
    with open(args.barcodes, "r") as f:
        barcodes = [line.strip() for line in f if line.strip()]

    # Building matrix

    sample_metrics = pl.read_parquet(args.parquet).to_pandas().set_index("CB").loc[barcodes]

    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = args.fragments,
        path_to_regions = args.regions,
        path_to_blacklist = args.blacklist,
        metrics = sample_metrics,
        valid_bc = barcodes,
        n_cpu = args.n_cpu,
        project = args.sample,
        split_pattern = '-'
    )
    
    print('Before filtering peaks:')
    print(cistopic_obj)
    print()

    ### Manually subsetting peaks
    
    # Create mask
    mask_peaks = (
                (cistopic_obj.region_data['cisTopic_nr_frag'] > args.min_fragments_per_region) &
                (cistopic_obj.region_data['cisTopic_nr_acc'] > args.min_cells_per_region)
            )
    # Filter region data (var) and fragment matrix (X)
    cistopic_obj.region_data = cistopic_obj.region_data.loc[mask_peaks]
    cistopic_obj.fragment_matrix = cistopic_obj.fragment_matrix[mask_peaks.values, :]
    cistopic_obj.binary_matrix = cistopic_obj.binary_matrix[mask_peaks.values, :]
    # update internal n_regions and region names
    cistopic_obj.n_regions = cistopic_obj.region_data.shape[0]
    cistopic_obj.region_names = cistopic_obj.region_data.index.tolist()

    print('Before filtering peaks:')
    print(cistopic_obj)
    print()

    # Write to file
    pickle.dump(
        cistopic_obj,
        open(args.out_cistopic_obj, "wb")
        )

if __name__ == "__main__":
    main()
