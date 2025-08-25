import os
import argparse
import timeit
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import pandas as pd

def write_data(adata,annotations,output_file):
    # Start with obs
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

    df.to_csv(annotations, sep=",", index=True)
    adata.write_h5ad(output_file)

def initialize_parser():
    parser = argparse.ArgumentParser(description='Correcting batches and samples using Harmony')
    parser.add_argument('--input', type=str, required=True, help="h5ad list")
    parser.add_argument('--output', type=str, required=True, help="h5ad")
    parser.add_argument('--annotations', type=str, required=True, help="csv")
    parser.add_argument('--plotdir', type=str, required=True)
    parser.add_argument('--celltype_key', type=str, required=True)
    parser.add_argument('--sample_key', type=str, required=True)
    parser.add_argument('--batch_key', type=str, required=True)
    return parser

# -----------------------------
# Main
# -----------------------------
def main():
    parser = initialize_parser()
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    annotations = args.annotations
    plotdir = args.plotdir
    celltype_key = args.celltype_key
    batch_key = args.batch_key
    sample_key = args.sample_key

    sc.settings.figdir = plotdir

    print("1. Load Data")
    start = timeit.default_timer()
    adata = ad.io.read_h5ad(input_file)
    stop = timeit.default_timer()
    print(f"Data loaded in {round(stop-start,2)}s")
    

    print("2. Batch correction (Harmony)")
    start = timeit.default_timer()
    sce.pp.harmony_integrate(adata, key=[batch_key,sample_key])
    stop = timeit.default_timer()
    print(f"Batches corrected in {round(stop-start,2)}s")

    print("4. Plotting UMAPs")
    start = timeit.default_timer()
    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    sc.tl.umap(adata)
    sc.pl.umap(adata,
                color=[sample_key,batch_key],
                show=False,
                save=f"_batch_corrected_merge.png")
    sc.pl.umap(adata,
                color=celltype_key,
                show=False,
                save=f"_annotation_batch_corrected_merge.png")
    stop = timeit.default_timer()
    print(f"UMAPs plotted in {round(stop-start,2)}s")

    print("5. Write data")
    start = timeit.default_timer()
    write_data(adata,annotations,output_file)
    stop = timeit.default_timer()
    print(f"Data written in {round(stop-start,2)}s")

    
if __name__ == "__main__":
    start_total = timeit.default_timer()
    main()
    stop_total = timeit.default_timer()
    print(f"Total runtime: {round(stop_total-start_total,2)}s")
