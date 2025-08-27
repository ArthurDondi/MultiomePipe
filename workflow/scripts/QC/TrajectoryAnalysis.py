import os
import json
import argparse
import timeit
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd



def initialize_parser():
    parser = argparse.ArgumentParser(description='Correcting batches and samples using Harmony')
    parser.add_argument('--input', type=str, required=True, help="h5ad")
    parser.add_argument('--output', type=str, required=True, help="h5ad")
    parser.add_argument('--diffusion_component', type=int, required=True)
    parser.add_argument('--plotdir', type=str, required=True)
    return parser

# -----------------------------
# Main
# -----------------------------
def main():
    parser = initialize_parser()
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    diffusion_component = args.diffusion_component 
    plotdir = args.plotdir

    sc.settings.figdir = plotdir

    print("1. Load Data")
    start = timeit.default_timer()
    adata = ad.io.read_h5ad(input_file)
    stop = timeit.default_timer()
    print(f"Data loaded in {round(stop-start,2)}s")
    

    print("2. Initializing root")
    start = timeit.default_timer()
    
    # Initializing root on user defined component
    root_ixs = adata.obsm["X_diffmap"][:, diffusion_component].argmin()
    adata.uns["iroot"] = root_ixs
    
    # Computing diffusion pseudotime
    sc.tl.dpt(adata)
    sc.pl.umap(adata, 
               color='dpt_pseudotime',
               show=False,
               save=f"_diffusion_pseudotime_merge.png")
    
    # Getting top dynamic genes
    pt = adata.obs['dpt_pseudotime']
    corrs = []
    for g in adata.var_names:
        corrs.append(np.corrcoef(adata[:, g].X.toarray().flatten(), pt)[0,1])

    adata.var['pseudotime_corr'] = corrs
    dynamic_genes = adata.var.sort_values('pseudotime_corr', ascending=False).head(500).index

    #sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata[:, dynamic_genes])

    # Compute gene-gene neighborhood graph
    sc.pp.neighbors(adata[:, dynamic_genes], use_rep='X_pca')
    sc.tl.leiden(adata, resolution=1.0, key_added='gene_modules')

    sc.pl.dpt_groups_pseudotime(adata,
                                show=False,
                                save=f"_merge.png"
                                )

    # Heatmap
    sc.pl.heatmap(
        adata,
        var_names=dynamic_genes,
        groupby='gene_modules',
        use_raw=False,
        swap_axes=True,
        show=False,
        save=f"_heatmap_merge.png")

    stop = timeit.default_timer()
    print(f"Batches corrected in {round(stop-start,2)}s")


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
