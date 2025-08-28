import os
import json
import argparse
import timeit
import anndata as ad
import scanpy as sc
import numpy as np
import scanpy as sc
import numpy as np
import pandas as pd
from tqdm import tqdm 
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


def compute_pseudotime(adata,celltype_key,root_ctype,diffusion_component):
    ### Initializing root on user defined component
    # subset of precursor cells
    precursors = adata[adata.obs[celltype_key]==root_ctype]

    # pick the cell with minimal value along the chosen diffusion component
    root_cell = precursors.obs_names[precursors.obsm["X_diffmap"][:, diffusion_component].argmin()]

    # get its integer position in the full adata
    adata.uns["iroot"] = list(adata.obs_names).index(root_cell)

    ### Computing diffusion pseudotime
    sc.tl.dpt(adata,n_branchings=1)
    # Computing PAGA connectivity graph
    sc.tl.paga(adata, groups='cell_type')
    # Plot PAGA first, so that `adata.uns['paga']['pos']` exists.
    sc.pl.paga(adata, 
               color=['cell_type', 'dpt_pseudotime'],
               show=False,
               save=f"_merge.png")
    # Recomputing UMAP with PAGA initialization
    sc.tl.umap(adata, init_pos='paga')
    
    return adata
    

def identify_dynamic_genes(adata,quantile):
    
    pseudotime = adata.obs['dpt_pseudotime'].values
    results = []

    adata_hvg = adata[:,adata.var['highly_variable']]
    
    for gene in tqdm(adata_hvg.var_names, desc="Processing genes"):
        expr = adata_hvg[:, gene].X.toarray().flatten()
        rho, pval = spearmanr(expr, pseudotime)
        results.append([gene, rho, pval])

    # Put into dataframe
    df = pd.DataFrame(results, columns=["gene", "rho", "pval"])

    # Bonferroni correction
    df["pval_adj"] = multipletests(df["pval"], method="bonferroni")[1]

    # Add -log10(pval_adj) for convenience
    df["neglog10_pval"] = -np.log10(df["pval_adj"].replace(0, np.nextafter(0,1)))

    # Rank genes by adjusted significance + effect size
    df["score"] = df["neglog10_pval"] * df["rho"].abs()

    df.sort_values("score", ascending=False, inplace=True)

    threshold = df["score"].quantile(quantile)  
    dynamic_genes = df[df["score"] >= threshold]['gene']

    adata_dyn = adata[:, dynamic_genes].copy()

    return adata_dyn,dynamic_genes

def plotting(adata,adata_dyn,dynamic_genes,celltype_key,cat_order):

    adata_dyn.obs['cell_type']= pd.Categorical(adata_dyn.obs['cell_type'],
                           categories=cat_order,
                           ordered=True )
    
    cell_order = adata_dyn.obs.sort_values(['cell_type', 'dpt_pseudotime']).index
    adata_dyn = adata_dyn[cell_order, :].copy()

    sc.pl.heatmap(
        adata_dyn,
        var_names=dynamic_genes,
        groupby='cell_type',
        use_raw=False,
        swap_axes=True,
        show=False,
        save=f"_cell_type_merge.png")

    sc.pl.heatmap(
        adata_dyn,
        var_names=dynamic_genes,
        groupby='dpt_pseudotime',
        use_raw=False,
        swap_axes=True,
        show=False,
        save=f"_dpt_pseudotime_merge.png")
    
    sc.pl.scatter(adata,
                  basis="diffmap",
                  color=celltype_key,
                  components=[2, 3],
                  show=False,
                  save=f"_diffusionmap_DC2-3_merge.png")
    
    sc.pl.scatter(adata,
                  basis="diffmap",
                  color=celltype_key,
                  components=[1, 2],
                  show=False,
                  save=f"_diffusionmap_DC1-2_merge.png")
    
    sc.pl.umap(adata, 
            color=['dpt_pseudotime','cell_type'],
            edges=True,
            show=False,
            save=f"_diffusion_pseudotime_merge.png")
    sc.pl.umap(adata, 
            color='dpt_groups',
            edges=True,
            show=False,
            save=f"_dpt_groups_merge.png")

def initialize_parser():
    parser = argparse.ArgumentParser(description='Pseudotime analysis using sc.tl.dpt')
    parser.add_argument('--input', type=str, required=True, help="h5ad")
    parser.add_argument('--output', type=str, required=True, help="h5ad")
    parser.add_argument('--diffusion_component', type=int, required=True)
    parser.add_argument('--quantile', type=float, required=True, default = 0.99)
    parser.add_argument('--celltype_key', type=str, required=True)
    parser.add_argument('--root_ctype', type=str, required=True)
    parser.add_argument('--plotdir', type=str, required=True)
    parser.add_argument('--celltype_mask', type=str,  nargs='+', required=True)
    parser.add_argument('--cat_order', type=str,  nargs='+', required=True)
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
    root_ctype = args.root_ctype
    celltype_key = args.celltype_key
    celltype_mask = args.celltype_mask
    quantile = args.quantile
    cat_order = args.cat_order

    sc.settings.figdir = plotdir

    print("1. Load Data")
    start = timeit.default_timer()
    adata = ad.io.read_h5ad(input_file)
    # Keep only cell types of interest:
    adata = adata[~adata.obs['cell_type'].isin(celltype_mask)]
    stop = timeit.default_timer()
    print(f"Data loaded in {round(stop-start,2)}s")
    

    print("2. Computing diffusion pseudotime")
    start = timeit.default_timer()
    adata = compute_pseudotime(adata,celltype_key,root_ctype,diffusion_component)
    stop = timeit.default_timer()
    print(f"Computing diffusion pseudotime in {round(stop-start,2)}s")

    print("3. Identifying dynamic genes")
    start = timeit.default_timer()
    adata_dyn, dynamic_genes = identify_dynamic_genes(adata,quantile)
    stop = timeit.default_timer()
    print(f"Dynamic genes identified in {round(stop-start,2)}s")

    print("4. Plotting")
    start = timeit.default_timer()
    plotting(adata,adata_dyn,dynamic_genes,celltype_key,cat_order)
    stop = timeit.default_timer()
    print(f"Dynamic genes identified in {round(stop-start,2)}s")
    
    
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
