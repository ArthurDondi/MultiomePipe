import os
import argparse
import timeit
import anndata as ad
import scanpy as sc
import json
import pandas as pd

def plot_and_add_annotations(adata, annotation_file, celltype_key, doublets, leiden_res, mode):
    
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

    adata = plot_and_add_annotations(adata, annotation_file, celltype_key, doublets, leiden_res, mode)

    adata = prepare_trajectory_analysis(adata,celltype_key)

    write_data(adata,annotation_csv,output_file)
    
if __name__ == "__main__":
    start_total = timeit.default_timer()
    main()
    stop_total = timeit.default_timer()
    print(f"Total runtime: {round(stop_total-start_total,2)}s")
