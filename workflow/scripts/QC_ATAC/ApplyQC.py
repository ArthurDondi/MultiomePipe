#!/usr/bin/env python3
r"""
Apply QC filters to a single sample and generate QC plots.

This script:
- Computes barcodes passing QC thresholds (user-defined or automatic)
- Generates sample-level statistics plots
- Generates barcode-level statistics plots (with detailed info)
- Writes the barcodes passing QC to a text file (one barcode per line)

Usage:
python ApplyQC.py --outdir outs/qc_output \
    --sample sampleA \
    --barcodes barcodes_passing_qc.tsv \
    [--unique_fragments_threshold 1000] \
    [--tss_enrichment_threshold 5] \
    [--frip_threshold 0.2]
"""

import argparse
from pycisTopic.plotting.qc_plot import plot_sample_stats, plot_barcode_stats
from pycisTopic.qc import get_barcodes_passing_qc_for_sample
import matplotlib
matplotlib.use("Agg") 


def none_or_float(x):
    if x in ("None", "none", ""):
        return None
    return float(x)

    
def none_or_int(x):
    if x in ("None", "none", ""):
        return None
    return int(x)



def main():
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", type = str, required=True)
    p.add_argument("--sample", type = str, required=True)
    p.add_argument("--barcodes", required = True)
    p.add_argument("--unique_fragments_threshold", type=none_or_int, default=None)
    p.add_argument("--tss_enrichment_threshold", type=none_or_float, default=None)
    p.add_argument("--frip_threshold", type = int, default=0)

    args = p.parse_args()

    plotdir = f"{args.outdir}/Plots"

    fig = plot_sample_stats(sample_id = args.sample, pycistopic_qc_output_dir = args.outdir)

    fig.savefig(f"{plotdir}/sample_statistics.png", dpi=300)

    barcodes_passing_filters, thresholds= get_barcodes_passing_qc_for_sample(
            sample_id = args.sample,
            pycistopic_qc_output_dir = args.outdir,
            unique_fragments_threshold = args.unique_fragments_threshold,
            tss_enrichment_threshold = args.tss_enrichment_threshold, 
            frip_threshold = args.frip_threshold,
            use_automatic_thresholds = True, # is ignored for unique_fragments_threshold or tss_enrichment_threshold when not None
    )

    fig = plot_barcode_stats(sample_id = args.sample,
                             pycistopic_qc_output_dir = args.outdir,
                             bc_passing_filters = barcodes_passing_filters,
                             detailed_title = False,
                             **thresholds)
    
    fig.savefig(f"{plotdir}/cells_statistics.png", dpi=300)

    fig = plot_barcode_stats(sample_id = args.sample,
                             pycistopic_qc_output_dir = args.outdir,
                             bc_passing_filters = barcodes_passing_filters,
                             detailed_title = True,
                             **thresholds)
    
    fig.savefig(f"{plotdir}/cells_statistics_info.png", dpi=300)

    with open(args.barcodes, "w") as f:
        for barcode in barcodes_passing_filters:
            f.write(barcode + "\n")

if __name__ == "__main__":
    main()
