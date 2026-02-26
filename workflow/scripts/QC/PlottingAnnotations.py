import os
import argparse
import timeit
import anndata as ad
import scanpy as sc
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_and_add_annotations(adata, annotation_file, plotdir, celltype_key, doublets, leiden_res, mode):    
    
    if mode == 'manual':
        with open(annotation_file, "r") as f:
            annotations = json.load(f)

        adata.obs[celltype_key] = adata.obs[leiden_res].map(annotations)

    sc.pl.umap(adata,
            color=celltype_key,
            legend_loc="on data",
            legend_fontsize=8,
            show=False,
            save=f"_annotated.png")
    
    adata = adata[~adata.obs[celltype_key].isin(doublets)].copy()

    sc.pl.umap(adata,
            color=celltype_key,
            legend_loc="on data",
            legend_fontsize=8,
            show=False,
            save=f"_annotated_clustersfiltered.png")
    
    # Removing ribosomal and mitochondrial genes for DEA
    mask_var = ~(
                adata.var_names.str.startswith(("MT-", "mt-")) |
                adata.var_names.str.startswith(("RPS", "RPL"))
                )

    sc.tl.rank_genes_groups(adata, groupby=celltype_key, mask_var = mask_var, method="wilcoxon")
    sc.pl.rank_genes_groups_dotplot(adata, 
                                    groupby=celltype_key,
                                    standard_scale="var",
                                    n_genes=5,
                                    show=False,
                                    swap_axes=True,
                                    save=f"_diff_markers_celltype_merge.png")
    
    ### Plot cell types distribution
    adata.obs[celltype_key] = adata.obs[celltype_key].astype("category")
    counts = adata.obs[celltype_key].value_counts().reindex(adata.obs[celltype_key].cat.categories)
    plt.bar(
        [0], 
        counts, 
        bottom=counts.cumsum().shift(fill_value=0), 
        color=adata.uns.get(celltype_key +"_colors", None)
    )
    plt.xticks([0], [celltype_key], rotation=45, ha='right')
    plt.ylabel("Number of cells")
    plt.legend(title="Cell type", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.title("Stacked histogram of cell types")
    plt.tight_layout()
    plt.savefig(f"{plotdir}/stacked_histogram_global.png", dpi=300, bbox_inches="tight")
    plt.close()

    ### Plot cell types distribution per sample
    ax = sns.histplot(
        data=adata.obs,
        x="sample",
        hue=celltype_key,
        multiple="stack",
        palette=adata.uns.get(f"{celltype_key}_colors", None),
        legend=False
    )

    plt.ylabel("Number of cells")
    plt.xticks(rotation=45, ha='right')
    plt.title("Stacked histogram of cell types per sample")

    # Manually create legend
    categories = adata.obs[celltype_key].cat.categories if hasattr(adata.obs[celltype_key], "cat") else adata.obs[celltype_key].unique()
    colors = adata.uns.get(f"{celltype_key}_colors", None)
    if colors is None:
        colors = sns.color_palette("tab20", len(categories)).as_hex()
    color_map = dict(zip(categories, colors))

    handles = [plt.Rectangle((0,0),1,1,color=color_map[cat]) for cat in categories]
    plt.legend(handles, categories, title="Cell type", bbox_to_anchor=(1.05,1), loc="upper left")

    plt.tight_layout()
    plt.savefig(f"{plotdir}/stacked_histogram_per_sample.png", dpi=300, bbox_inches="tight")
    plt.close()
    
    return adata

def prepare_trajectory_analysis(adata,celltype_key):
    sc.tl.diffmap(adata)
    sc.pl.scatter(adata,
                  basis="diffmap",
                  color=celltype_key,
                  components=[2, 3],
                  show=False,
                  save=f"_diffusionmap2-3_merge.png")
    sc.pl.scatter(adata,
                  basis="diffmap",
                  color=celltype_key,
                  components=[1, 2],
                  show=False,
                  save=f"_diffusionmap1-2_merge.png")

    return adata

def write_data(adata,annotation_csv,output_file):
    
    df = adata.obs.copy()

    # Add one embedding (e.g. UMAP)
    if "X_umap" in adata.obsm:
        umap_df = pd.DataFrame(
            adata.obsm["X_umap"],
            index=adata.obs.index,
            columns=["UMAP1", "UMAP2"]
        )
        df = pd.concat([df, umap_df], axis=1)

    if "X_umap_nocorrection" in adata.obsm:
        umap_df = pd.DataFrame(
            adata.obsm["X_umap_nocorrection"],
            index=adata.obs.index,
            columns=["UMAP1_nocorrection", "UMAP2_nocorrection"]
        )
        df = pd.concat([df, umap_df], axis=1)

    df.to_csv(annotation_csv, sep=",", index=True)
    adata.write_h5ad(output_file)

def initialize_parser():
    parser = argparse.ArgumentParser(description='QC + clustering for 10X Multiome RNA')
    parser.add_argument('--input', type=str, required=True, help="h5ad")
    parser.add_argument('--output', type=str, required=True, help="h5ad")
    parser.add_argument('--annotation_input', type=str, default='', help="json")
    parser.add_argument('--annotation_output', type=str, help="csv")
    parser.add_argument('--plotdir', type=str, required=True)
    parser.add_argument('--doublets', type=str, nargs='+', required=True, help="List of names of clusters to remove")
    parser.add_argument('--mode', type=str, required=True, choices=['auto', 'manual'], help="Automatic or manual annotation")
    parser.add_argument('--celltype_key', type=str, default ="cell_type", help="Name of celltype .obs column, default: 'cell_type'")
    parser.add_argument('--leiden_res', type=str, default="leiden_res_1.00", help="leiden resolution used for manual annotation")
    return parser

# -----------------------------
# Main
# -----------------------------
def main():
    parser = initialize_parser()
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    annotation_file = args.annotation_input
    annotation_csv = args.annotation_output
    plotdir = args.plotdir
    doublets = args.doublets
    mode = args.mode
    celltype_key = args.celltype_key
    leiden_res = args.leiden_res

    sc.settings.figdir = plotdir

    adata = ad.io.read_h5ad(input_file)

    adata = plot_and_add_annotations(adata, annotation_file, plotdir, celltype_key, doublets, leiden_res, mode)

    adata = prepare_trajectory_analysis(adata,celltype_key)

    write_data(adata,annotation_csv,output_file)
    
if __name__ == "__main__":
    start_total = timeit.default_timer()
    main()
    stop_total = timeit.default_timer()
    print(f"Total runtime: {round(stop_total-start_total,2)}s")
