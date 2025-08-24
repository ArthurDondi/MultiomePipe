import os
import argparse
import timeit
import anndata as ad
import scanpy as sc
import json

# -----------------------------
# Functions
# -----------------------------
def merge_data(input_files):
    adatas = [ad.io.read_h5ad(f) for f in input_files]
    for adata in adatas:
        adata.var_names_make_unique()
    print("Concat")
    adata_concat = ad.concat(
        adatas,
        merge="first",
        uns_merge="first",
        join="outer",
        label="split",
        keys=[f"split{i}" for i in range(len(adatas))]
    )
    print("Concat done")
    del adatas
    adata_concat.X = adata_concat.layers['counts']
    return adata_concat

def plot_qc_metrics(adata):

    print("Plotting violin")
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
        save=f"_merge.png"
    )

    print("Plotting scatter")
    sc.pl.scatter(
        adata,
        "total_counts",
        "n_genes_by_counts",
        color="pct_counts_mt",
        size=50,
        show=False,
        save=f"_merge.png"
    )


def run_normalization_and_clustering(adata, celltype_key, batch_key, sample_key, is_filtered):

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    # Plot highly variable genes:
    sc.pl.highly_variable_genes(adata,
                                show=False,
                                save=f"_HVG_merge.png")

    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    # Plot PCs variance ratio :
    sc.pl.pca_variance_ratio(adata, 
                             n_pcs=50, 
                             log=True,
                             show=False,
                             save=f"_merge.png")
    sc.tl.umap(adata)
    for res in [0.2, 0.5, 1]:
        sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res) # cluster labels in adata.obs['leiden']
    sc.pl.umap(adata,
            color=[f"leiden_res_{res:4.2f}" for res in [0.2, 0.5, 1]],
            legend_loc="on data",
            show=False,
            save=f"_leiden_res_merge.png")
    
    if is_filtered:
        QCs = ["log1p_total_counts", "pct_counts_mt", "doublet_score"]
    else:
        QCs = ["log1p_total_counts", "pct_counts_mt", "doublet_score", "background_fraction"]
    sc.pl.umap(adata,
                color=QCs,
                show=False,
                save=f"_QC_merge.png")
    sc.pl.umap(adata,
                color=[batch_key,sample_key],
                show=False,
                save=f"_batch_nocorrection_merge.png")
    sc.pl.umap(adata,
                color=celltype_key,
                legend_loc="on data",
                legend_fontsize=8,
                show=False,
                save=f"_annotated_merge.png")
    
    return adata

def prepare_annotation(adata,markers_file):

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

    for res in [0.2, 0.5, 1]:
        sc.pl.dotplot(
            adata,
            groupby=f"leiden_res_{res:4.2f}",
            var_names=marker_genes_in_data,
            standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
            show=False,
            save=f"_leiden_res_{res:4.2f}_merge.png"
        )
    
    return adata

def initialize_parser():
    parser = argparse.ArgumentParser(description='Merging adatas from all samples and plotting')
    parser.add_argument('--input', type=str, nargs='+', required=True, help="h5ad list")
    parser.add_argument('--output', type=str, required=True, help="h5ad")
    parser.add_argument('--markers', type=str, required=True, help="json")
    parser.add_argument('--plotdir', type=str, required=True)
    parser.add_argument('--celltype_key', type=str, required=True)
    parser.add_argument('--sample_key', type=str, required=True)
    parser.add_argument('--batch_key', type=str, required=True)
    parser.add_argument('--is_filtered', action='store_true')
    return parser

# -----------------------------
# Main
# -----------------------------
def main():
    parser = initialize_parser()
    args = parser.parse_args()

    input_files = args.input
    output_file = args.output
    markers = args.markers
    plotdir = args.plotdir
    celltype_key = args.celltype_key
    batch_key = args.batch_key
    sample_key = args.sample_key
    is_filtered = args.is_filtered

    sc.settings.figdir = plotdir

    # 1. Merge data
    print("1. Merge data")
    start = timeit.default_timer()
    adata = merge_data(input_files)
    stop = timeit.default_timer()
    print(f"Loaded data in {round(stop-start,2)}s")

    # 2. QC metrics
    print("2. QC metrics")
    start = timeit.default_timer()
    plot_qc_metrics(adata)
    stop = timeit.default_timer()
    print(f"QC metrics computed in {round(stop-start,2)}s")

    # 3. Normalization + clustering
    print("3. Normalization + clustering")
    start = timeit.default_timer()
    adata = run_normalization_and_clustering(adata, celltype_key, batch_key, sample_key, is_filtered)
    stop = timeit.default_timer()
    print(f"Normalization + clustering done in {round(stop-start,2)}s")

    # 4. Annotation (run only if a file path is provided)
    print("4. Annotation")
    if markers != '':
        start = timeit.default_timer()
        adata = prepare_annotation(adata,markers)
        stop = timeit.default_timer()
        print(f"Annotation done in {round(stop-start,2)}s")
    else:
        print("No marker genes file provided")

    # 5. Write data
    print("5. Write data")
    start = timeit.default_timer()
    adata.write_h5ad(output_file)
    stop = timeit.default_timer()
    print(f"Data written in {round(stop-start,2)}s")

if __name__ == "__main__":
    start_total = timeit.default_timer()
    main()
    stop_total = timeit.default_timer()
    print(f"Total runtime: {round(stop_total-start_total,2)}s")
