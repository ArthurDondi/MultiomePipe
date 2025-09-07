#!/usr/bin/env python3
r"""
Apply QC filters to a single sample and generate QC plots.

This script:
- Loads a pickled cisTopic object
- Runs LDA models with user-defined parameters
- Evaluates models and saves evaluation plot
- Writes updated cisTopic object back to disk

Usage:
python CreateATACCountMatrix.py \
    --outdir outs/qc_output \
    --in_cistopic_obj input.pkl \
    --out_cistopic_obj output.pkl \
    --n_cpu 8 \
    --n_topics 2 5 10 15 20 \
    --n_iter 500 \
    --random_state 555 \
    --alpha 50 \
    --alpha_by_topic \
    --eta 0.1 \
    --eta_by_topic
"""

import argparse
from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg") 
import pickle


def main():
    p = argparse.ArgumentParser(description="Apply QC filters and build cisTopic object.")
    p.add_argument("--outdir", type=str, required=True, help="Output directory for QC results.")
    p.add_argument("--in_cistopic_obj", required=True, help="Path to input pickled cistopic object.")
    p.add_argument("--out_cistopic_obj", required=True, help="Path to output pickled cistopic object.")
    p.add_argument("--cell_data", required=True, help="cell_data.tsv (index is barcode)")
    p.add_argument("--sample", required=True)
    p.add_argument("--sample_key", default="sample", help="sample id column in cell_data")

    args = p.parse_args()

    plotdir = f"{args.outdir}/Plots"

    # Load input cisTopic object
    with open(args.in_cistopic_obj, "rb") as f:
        cistopic_obj = pickle.load(f)

    # read cell metadata
    cell_data = pd.read_csv(args.cell_data, index_col=0)
    cell_data = cell_data[cell_data[args.sample_key]==args.sample]
    # chabge BC-1_Sample by BC-1-Sample to match cistopic_obj.cell_data.index
    cell_data.index = cell_data.index.str.replace('-1_', '-1-')
    cell_data = cell_data.reindex(cistopic_obj.cell_data.index)

    cistopic_obj.add_cell_data(cell_data, split_pattern='-')

    cistopic_obj.cell_data = cistopic_obj.cell_data.fillna("NA")

    find_clusters(
        cistopic_obj,
        target  = 'cell',
        k = 10,
        res = [0.6, 1.2, 3],
        prefix = 'pycisTopic_',
        scale = True,
        split_pattern = '-'
    )

    run_umap(
    cistopic_obj,
    target  = 'cell', scale=True)

    run_tsne(
    cistopic_obj,
    target  = 'cell', scale=True)

    plot_metadata(
        cistopic_obj,
        reduction_name='UMAP',
        variables=['cell_type', 'pycisTopic_leiden_10_0.6', 'pycisTopic_leiden_10_1.2', 'pycisTopic_leiden_10_3'], 
        target='cell', num_columns=4,
        text_size=20,
        dot_size=5,
        save = f"{plotdir}/umap_ctypes_leiden.png")
    
    plot_metadata(
        cistopic_obj,
        reduction_name='UMAP',
        variables=['log10_unique_fragments_count', 'tss_enrichment', 'fraction_of_fragments_in_peaks'],
        target='cell', num_columns=4,
        text_size=20,
        dot_size=5,
        save = f"{plotdir}/umap_QC.png")
    
    plot_topic(
        cistopic_obj,
        reduction_name = 'UMAP',
        target = 'cell',
        num_columns=5, 
        save = f"{plotdir}/umap_topics.png")
    
    cell_topic_heatmap(
        cistopic_obj,
        variables = ['cell_type'],
        scale = False,
        legend_loc_x = 1.0,
        legend_loc_y = -1.2,
        legend_dist_y = -1,
        figsize = (10, 10),
        save = f"{plotdir}/heatmap.png"
    )

    # Write updated object
    with open(args.out_cistopic_obj, "wb") as f:
        pickle.dump(cistopic_obj, f)


if __name__ == "__main__":
    main()
