import os
import argparse
import timeit
import anndata as ad
import scanpy as sc
import json

# -----------------------------
# Functions
# -----------------------------
def load_data(input_file):
    adata = ad.io.read_h5ad(input_file)
    adata.var_names_make_unique()
    return adata 

def run_calculate_qc_metrics(adata, sample):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","ribo","hb"], inplace=True, log1p=True)
    adata.var_names_make_unique()

    # Plotting
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
        save=f"_{sample}.png"
    )
    sc.pl.scatter(
        adata,
        "total_counts",
        "n_genes_by_counts",
        color="pct_counts_mt",
        size=50,
        show=False,
        save=f"_{sample}.png"
    )
    return adata

def run_raw_filtering(adata, min_genes, min_cells, sample):
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Plotting# Plotting
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
        save=f"_{min_genes}minGenes_{min_cells}minCells_{sample}.png"
    )
    sc.pl.scatter(
        adata,
        "total_counts",
        "n_genes_by_counts",
        color="pct_counts_mt",
        size=50,
        show=False,
        save=f"_{min_genes}minGenes_{min_cells}minCells_{sample}.png"
    )
    return adata

def identify_doublets(adata):
    sc.pp.scrublet(adata)
    n_predicted_doublets = adata.obs['predicted_doublet'].sum()
    total_cells = adata.n_obs
    print(f"Predicted doublets: {n_predicted_doublets} over {total_cells} cells ({100 * n_predicted_doublets / total_cells:.1f}%)")
    return adata


def run_normalization_and_clustering(adata, sample):
    # Saving count data
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    # Plot highly variable genes:
    sc.pl.highly_variable_genes(adata,
                                show=False,
                                save=f"_HVG_{sample}.png")

    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    # Plot PCs variance ratio :
    sc.pl.pca_variance_ratio(adata, 
                             n_pcs=50, 
                             log=True,
                             show=False,
                             save=f"_{sample}.png")
    sc.tl.umap(adata)
    for res in [0.2, 0.5, 1]:
        sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res) # cluster labels in adata.obs['leiden']
    sc.pl.umap(adata,
            color=[f"leiden_res_{res:4.2f}" for res in [0.2, 0.5, 1]],
            show=False,
            save=f"_leiden_res_{sample}.png")
    sc.pl.umap(adata,
                color=["total_counts", "pct_counts_mt", "doublet_score", "background_fraction"],
                show=False,
                save=f"_QC_{sample}.png")

    return adata

def run_annotation(adata,markers_file,sample):

    with open(markers_file, "r") as f:
        marker_genes = json.load(f)

    # Keeping only genes in dataset
    marker_genes_in_data = {}
    for ct, markers in marker_genes.items():
        markers_found = []
        for marker in markers:
            if marker in adata.var.index:
                markers_found.append(marker)
        marker_genes_in_data[ct] = markers_found

    res = 0.5
    sc.pl.dotplot(
        adata,
        groupby=f"leiden_res_{res:4.2f}",
        var_names=marker_genes_in_data,
        standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
        show=False,
        save=f"_leiden_res_{res:4.2f}_{sample}.png"
    )
    
    return adata

def initialize_parser():
    parser = argparse.ArgumentParser(description='QC + clustering for 10X Multiome RNA')
    parser.add_argument('--input', type=str, required=True, help="h5ad")
    parser.add_argument('--output', type=str, required=True, help="h5ad")
    parser.add_argument('--markers', type=str, required=True, help="json")
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
    markers = args.markers
    sample = args.sample
    outdir = args.outdir
    min_genes = args.min_genes
    min_cells = args.min_cells

    sc.settings.figdir = outdir  + "/Plots"

    # 1. Load data
    print("1. Load data")
    start = timeit.default_timer()
    adata = load_data(input_file)
    stop = timeit.default_timer()
    print(f"Loaded data in {round(stop-start,2)}s")

    # 2. QC metrics
    print("2. QC metrics")
    start = timeit.default_timer()
    adata = run_calculate_qc_metrics(adata, sample)
    stop = timeit.default_timer()
    print(f"QC metrics computed in {round(stop-start,2)}s")

    # 3. RaW Filtering
    print("3. Filtering")
    start = timeit.default_timer()
    adata = run_raw_filtering(adata, min_genes, min_cells, sample)
    stop = timeit.default_timer()
    print(f"Filtering done in {round(stop-start,2)}s")

    # 4. Detecting Doublets
    print("4. Detecting Doublets")
    start = timeit.default_timer()
    adata = identify_doublets(adata)
    stop = timeit.default_timer()
    print(f"Detecting doublets done in {round(stop-start,2)}s")
    
    # 4. Normalization + clustering
    print("4. Normalization + clustering")
    start = timeit.default_timer()
    adata = run_normalization_and_clustering(adata, sample)
    stop = timeit.default_timer()
    print(f"Normalization + clustering done in {round(stop-start,2)}s")

    # 5. Annotation
    print("4. Annotation")
    start = timeit.default_timer()
    adata = run_annotation(adata,markers,sample)
    stop = timeit.default_timer()
    print(f"Annotation done in {round(stop-start,2)}s")

    # 6. Write data
    print("5. Write data")
    start = timeit.default_timer()
    adata.write(output_file)
    stop = timeit.default_timer()
    print(f"Data written in {round(stop-start,2)}s")

if __name__ == "__main__":
    start_total = timeit.default_timer()
    main()
    stop_total = timeit.default_timer()
    print(f"Total runtime: {round(stop_total-start_total,2)}s")
