import os
import argparse
import timeit
import gzip
import numpy as np
import pandas as pd
from scipy import io
import anndata as ad
import scanpy as sc
import muon as mu

# -----------------------------
# Functions
# -----------------------------
def load_data(input_file):
    mdata = mu.read_10x_h5(input_file)
    mdata.var_names_make_unique()
    if "rna" not in mdata.mod:
        raise ValueError(f"No 'rna' modality found in {input_file}. Available: {list(mdata.mod.keys())}")
    rna = mdata.mod['rna']
    return rna, mdata

def run_calculate_qc_metrics(rna, sample, outdir):
    rna.var["mt"] = rna.var_names.str.startswith("MT-")
    rna.var["ribo"] = rna.var_names.str.startswith(("RPS", "RPL"))
    rna.var["hb"] = rna.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(rna, qc_vars=["mt","ribo","hb"], inplace=True, log1p=True)
    rna.var_names_make_unique()

    # Plotting
    sc.pl.violin(
        rna,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
        save=f"{sample}_raw_violin.png"
    )
    sc.pl.scatter(
        rna,
        "total_counts",
        "n_genes_by_counts",
        color="pct_counts_mt",
        size=50,
        show=False,
        save=f"{sample}_raw_scatter.png"
    )
    return rna

def run_raw_filtering(rna, min_genes, min_cells, sample, outdir):
    sc.pp.filter_cells(rna, min_genes=min_genes)
    sc.pp.filter_genes(rna, min_cells=min_cells)

    # Plotting# Plotting
    sc.pl.violin(
        rna,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
        save=f"{sample}_{min_genes}minGenes_{min_cells}minCells_violin.png"
    )
    sc.pl.scatter(
        rna,
        "total_counts",
        "n_genes_by_counts",
        color="pct_counts_mt",
        size=50,
        show=False,
        save=f"{sample}_{min_genes}minGenes_{min_cells}minCells_scatter.png"
    )
    return rna

def run_normalization_and_clustering(rna, sample, outdir):
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna, n_top_genes=2000, subset=True)
    sc.pp.pca(rna, n_comps=30)
    sc.pp.neighbors(rna, n_neighbors=15, n_pcs=30)
    sc.tl.umap(rna)
    sc.tl.leiden(rna, key_added="clusters") # cluster labels in rna.obs['leiden']

    # Save clusters for SoupX
    cluster_file = os.path.join(outdir, f"{sample}_scanpy_clusters.tsv")
    rna.obs['clusters'].to_csv(cluster_file, sep='\t', header=False)
    print(f"Scanpy clusters saved to {cluster_file}")

    # Save cluster plots
    sc.pl.umap(rna, 
               color="clusters", 
               show=False,
               save=f"_{sample}_clusters.png")

    return rna

def write_10x_mtx(dirname, adata, var_names="gene_symbols", overwrite=False):
    if os.path.exists(dirname):
        if overwrite:
            import shutil
            shutil.rmtree(dirname)
        else:
            raise FileExistsError(f"{dirname} exists. Set overwrite=True to overwrite.")
    os.makedirs(dirname, exist_ok=True)

    io.mmwrite(os.path.join(dirname, "matrix.mtx"), adata.X.T)
    adata.obs.index.to_series().to_csv(
        os.path.join(dirname, "barcodes.tsv"), index=False, header=False
    )
    features = adata.var["gene_symbols"] if var_names=="gene_symbols" else adata.var[var_names]
    pd.Series(features).to_csv(
        os.path.join(dirname, "features.tsv"), index=False, header=False
    )

    for fname in ["matrix.mtx", "barcodes.tsv", "features.tsv"]:
        with open(os.path.join(dirname, fname), "rb") as f_in, gzip.open(os.path.join(dirname, fname+".gz"), "wb") as f_out:
            f_out.writelines(f_in)
        os.remove(os.path.join(dirname, fname))
    print(f"10X-style files written to {dirname}")

def write_data(rna, mdata, output, outdir, sample):
    # Muon H5
    mdata.write(output)
    # 10X export
    write_10x_mtx(
        os.path.join(outdir, "10x_for_SoupX"),
        rna,
        var_names="gene_ids",
        overwrite=True
    )

def initialize_parser():
    parser = argparse.ArgumentParser(description='QC + clustering for 10X Multiome RNA')
    parser.add_argument('--input', type=str, required=True)
    parser.add_argument('--output', type=str, required=True)
    parser.add_argument('--sample', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    parser.add_argument('--min_genes', type=int, default=100)
    parser.add_argument('--min_cells', type=int, default=3)
    return parser

# -----------------------------
# Main
# -----------------------------
def main():
    parser = initialize_parser()
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    sample = args.sample
    outdir = args.outdir
    min_genes = args.min_genes
    min_cells = args.min_cells

    os.makedirs(outdir, exist_ok=True)
    sc.settings.figdir = outdir  + "/Plots"

    # 1. Load data
    print("1. Load data")
    start = timeit.default_timer()
    rna, mdata = load_data(input_file)
    stop = timeit.default_timer()
    print(f"Loaded data in {round(stop-start,2)}s")

    # 2. QC metrics
    print("2. QC metrics")
    start = timeit.default_timer()
    rna = run_calculate_qc_metrics(rna, sample, outdir)
    stop = timeit.default_timer()
    print(f"QC metrics computed in {round(stop-start,2)}s")

    # 3. Filtering
    print("3. Filtering")
    start = timeit.default_timer()
    rna = run_raw_filtering(rna, min_genes, min_cells, sample, outdir)
    stop = timeit.default_timer()
    print(f"Filtering done in {round(stop-start,2)}s")

    # 4. Normalization + clustering
    print("4. Normalization + clustering")
    start = timeit.default_timer()
    rna = run_normalization_and_clustering(rna, sample, outdir)
    stop = timeit.default_timer()
    print(f"Normalization + clustering done in {round(stop-start,2)}s")

    # 5. Write data
    print("5. Write data")
    start = timeit.default_timer()
    write_data(rna, mdata, output_file, outdir, sample)
    stop = timeit.default_timer()
    print(f"Data written in {round(stop-start,2)}s")

if __name__ == "__main__":
    start_total = timeit.default_timer()
    main()
    stop_total = timeit.default_timer()
    print(f"Total runtime: {round(stop_total-start_total,2)}s")
