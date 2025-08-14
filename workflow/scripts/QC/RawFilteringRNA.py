import argparse
import timeit

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import muon as mu
import pysam

def run_calculate_qc_metrics(rna, sample):
    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    rna.var["mt"] = rna.var_names.str.startswith("MT-")
    # ribosomal genes
    rna.var["ribo"] = rna.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    rna.var["hb"] = rna.var_names.str.contains("^HB[^(P)]")
    # Calculate qc metrics
    sc.pp.calculate_qc_metrics(
        rna, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )
    rna.var_names_make_unique()
    # Plotting
    sc.pl.violin(
        rna,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
        save=f"_raw_{sample}.png"
    )
    sc.pl.scatter(rna, 
        "total_counts",
        "n_genes_by_counts",
        color="pct_counts_mt",
        size=50,
        show=False,
        save=f"_raw_{sample}.png"
    )
    return rna

def run_raw_filtering(rna, min_genes, min_cells, sample):
    # Filtering
    sc.pp.filter_cells(rna, min_genes=min_genes)
    sc.pp.filter_genes(rna, min_cells=min_cells)    
    # Plotting
    sc.pl.violin(
        rna,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
        save=f"{sample}_{min_genes}_min_genes_{min_cells}_min_cells_violin.png"
    )
    sc.pl.scatter(rna, 
        "total_counts",
        "n_genes_by_counts",
        color="pct_counts_mt",
        size=50,
        show=False,
        save=f"{sample}_{min_genes}_min_genes_{min_cells}_min_cells_pct_counts_mt.png"
    )
    return rna

def initialize_parser():
    parser = argparse.ArgumentParser(description='Script to perform initial QC on 10X Multiome GEX data')
    parser.add_argument('--input', type=str, help='Input file prefix', required = True)   
    parser.add_argument('--output', type=str, help='Output file prefix', required = True)
    parser.add_argument('--sample', type=str, help='Sample name', required = True)
    parser.add_argument('--plotdir', type=str,  help='Plotting directory', required = True)
    parser.add_argument('--min_genes', type=int, default = 100, help='Delta VAF between cancer and non-cancer cells', required = True)
    parser.add_argument('--min_cells', type=int, default = 3, help='Delta MCF (cancer cell fraction) between cancer and non-cancer cells', required = True)
    return (parser)

def main():

    # 1. Arguments
    parser = initialize_parser()
    args = parser.parse_args()

    input = args.input
    output = args.output
    sample = args.sample
    plotdir = args.plotdir
    min_genes = args.min_genes
    min_cells = args.min_cells

    sc.settings.figdir = plotdir

    # 1.1: Load data
    print ('\n- Loading Data\n')
    start = timeit.default_timer()
    # Load RNA+ATAC matrices
    mdata = mu.read_10x_h5(input)
    mdata.var_names_make_unique()
    # Keep RNA only
    rna = mdata.mod['rna']
    stop = timeit.default_timer()
    print ('\nLoading Data computing time: ' + str(round(stop - start,2)) + ' seconds')

    # 1.2: Calculate QC metrics
    print ('\n- Calculate QC metrics \n')
    start = timeit.default_timer()
    rna = run_calculate_qc_metrics(rna, sample)
    stop = timeit.default_timer()
    print ('\nCalculate QC metrics computing time: ' + str(round(stop - start,2)) + ' seconds')

    # 1.3: Filter cells and genes
    print ('\n- Filter cells and genes \n')
    start = timeit.default_timer()
    rna = run_raw_filtering(rna, min_genes, min_cells, sample)
    stop = timeit.default_timer()
    print ('\nFilter cells and genes computing time: ' + str(round(stop - start,2)) + ' seconds')

    # 1.4 Writing to disk
    print ('\n- Writing to disk \n')
    start = timeit.default_timer()
    mdata.write(output)
    stop = timeit.default_timer()
    print ('\nWriting to disk computing time: ' + str(round(stop - start,2)) + ' seconds')


#-------------------------
# Running scRNA somatic variant calling
#-------------------------
if __name__ == '__main__':
    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    print ('\nTotal computing time: ' + str(round(stop - start,2)) + ' seconds')