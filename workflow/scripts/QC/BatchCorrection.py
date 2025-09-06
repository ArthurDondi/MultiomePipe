import os
import json
import argparse
import timeit
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import pandas as pd

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
            save=f"_manual_markers_leiden_res_{res:4.2f}_merge.png"
        )
    
    return adata

def clustering(adata,sample_key,donor_key,is_filtered):

    if "X_pca_harmony" in adata.obsm:
        sc.pp.neighbors(adata, use_rep='X_pca_harmony')
        sc.tl.umap(adata)

    for res in [0.2, 0.5, 1]:
        sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res)
        sc.tl.rank_genes_groups(adata, groupby=f"leiden_res_{res:4.2f}", method="wilcoxon")
        sc.pl.rank_genes_groups_dotplot(adata, 
                                        groupby=f"leiden_res_{res:4.2f}",
                                        standard_scale="var",
                                        n_genes=5,
                                        show=False,
                                        save=f"_diff_markers_leiden_res_{res:4.2f}_merge.png")
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
                color=[donor_key,sample_key],
                show=False,
                save=f"_merge.png")
    
    
    sc.pl.umap(adata,
            color=[f"leiden_res_{res:4.2f}" for res in [0.2, 0.5, 1]],
            legend_loc="on data",
            show=False,
            save=f"_leiden_res_merge.png")

def initialize_parser():
    parser = argparse.ArgumentParser(description='Correcting batches and samples using Harmony')
    parser.add_argument('--input', type=str, required=True, help="h5ad")
    parser.add_argument('--output', type=str, required=True, help="h5ad")
    parser.add_argument('--markers', type=str, required=True, help="json")
    parser.add_argument('--plotdir', type=str, required=True)
    parser.add_argument('--is_filtered', action='store_true')
    parser.add_argument('--sample_key', type=str, required=True)
    parser.add_argument('--donor_key', type=str, required=True)
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
    plotdir = args.plotdir
    donor_key = args.donor_key
    sample_key = args.sample_key
    is_filtered = args.is_filtered

    sc.settings.figdir = plotdir

    print("1. Load Data")
    start = timeit.default_timer()
    adata = ad.io.read_h5ad(input_file)
    stop = timeit.default_timer()
    print(f"Data loaded in {round(stop-start,2)}s")
    

    print("2. Batch correction (Harmony)")
    start = timeit.default_timer()
    if (len(adata.obs[donor_key].unique()) > 1) | (len(adata.obs[sample_key].unique()) > 1):
        sce.pp.harmony_integrate(adata, 
                                key=[donor_key,sample_key],
                                theta=3.0,           # default ~2; smaller = weaker integration
                                sigma=0.1)
    stop = timeit.default_timer()
    print(f"Batches corrected in {round(stop-start,2)}s")

    print("3. Clustering")
    start = timeit.default_timer()
    clustering(adata,sample_key,donor_key,is_filtered)
    stop = timeit.default_timer()
    print(f"Clustered in {round(stop-start,2)}s")

    print("4. Plot dotplots for Annotation")
    if markers != '':
        start = timeit.default_timer()
        adata = prepare_annotation(adata,markers)
        stop = timeit.default_timer()
        print(f"Dotplots plotted in in {round(stop-start,2)}s")
    else:
        print("No marker genes file provided")

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
