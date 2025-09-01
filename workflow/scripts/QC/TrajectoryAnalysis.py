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
from scipy.cluster.hierarchy import linkage, fcluster


def compute_pseudotime(adata,celltype_key,root_ctype,diffusion_component,branching):
    ### Initializing root on user defined component
    # subset of precursor cells
    precursors = adata[adata.obs[celltype_key]==root_ctype]

    # pick the cell with minimal value along the chosen diffusion component
    if len(precursors)>0:
        root_cell = precursors.obs_names[precursors.obsm["X_diffmap"][:, diffusion_component].argmin()]
    else:
        print("No precursor found, selecting a random cell")
        root_cell = np.random.choice(adata.obs_names)


    # get its integer position in the full adata
    adata.uns["iroot"] = list(adata.obs_names).index(root_cell)

    ### Computing diffusion pseudotime
    sc.tl.dpt(adata,n_branchings=branching)

    try:
        # Computing PAGA connectivity graph
        sc.tl.paga(adata, groups='cell_type')
        # Plot PAGA first, so that `adata.uns['paga']['pos']` exists.
        sc.pl.paga(adata, 
                color='dpt_pseudotime',
                show=False,
                save=f"_merge.png")
        # Recomputing UMAP with PAGA initialization

        sc.tl.umap(adata, init_pos='paga')
        return adata
    # Sometimes a bug occurs if the matrix is too sparse
    except AttributeError:
        return adata

def cluster_genes(adata,gene_list,clustering_distance):
    # Extract expression (cells x genes)
    expr = adata[:, gene_list].X
    if not isinstance(expr, np.ndarray):
        expr = expr.toarray()  # convert sparse to dense if needed

    # Convert to DataFrame for convenience
    df = pd.DataFrame(expr, columns=gene_list)

    # Compute gene-gene correlation matrix
    corr_matrix = df.corr()

    # Hierarchical clustering (distance = 1 - correlation)
    Z = linkage(1 - corr_matrix, method='average')

    # Assign clusters (e.g., 3 clusters)
    #n_clusters = len(adata_dyn.obs['cell_type'].cat.categories)
    cluster_labels = fcluster(Z, t=clustering_distance, criterion='distance')

    # Create dictionary: cluster -> list of genes
    gene_clusters = {}
    for label in np.unique(cluster_labels):
        genes_in_cluster = [gene for gene, l in zip(gene_list, cluster_labels) if l == label]
        gene_clusters[f'cluster_{label}'] = genes_in_cluster
    
    return gene_clusters

def identify_dynamic_genes(adata,n_genes,clustering_distance,cat_order,celltype_mask):
    
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

    dynamic_genes = df['gene'][:n_genes]

    adata_dyn = adata[:, dynamic_genes].copy()

    gene_clusters = cluster_genes(adata_dyn,dynamic_genes,clustering_distance)

    # Saving celltype to colors mapping before sorting
    old_cats = adata.obs['cell_type'].cat.categories
    old_colors = adata.uns['cell_type_colors']
    cat_to_color = dict(zip(old_cats, old_colors))
    
    #Sorting cells for nicer heatmap
    cat_order = [i  for i in cat_order if i not in celltype_mask]
    adata_dyn.obs['cell_type']= pd.Categorical(adata_dyn.obs['cell_type'],
                           categories=cat_order,
                           ordered=True )
    cell_order = adata_dyn.obs.sort_values(['cell_type', 'dpt_pseudotime']).index
    adata_dyn = adata_dyn[cell_order, :].copy()

    # Build color array in the new order to keep original cell UMAP coloring (used in heatmap)
    adata_dyn.uns['cell_type_colors'] = np.array([cat_to_color[ct] for ct in cat_order if ct in cat_to_color])

    return adata_dyn, gene_clusters

def plotting(adata,adata_dyn,gene_clusters,celltype_key,branching):

    sc.pl.heatmap(
        adata_dyn,
        var_names=gene_clusters,
        groupby='cell_type',
        use_raw=False,
        swap_axes=True,
        show=False,
        save=f"_dpt_ctype.png")
    
    sc.pl.heatmap(
        adata_dyn,
        var_names=gene_clusters,
        groupby='dpt_pseudotime',
        use_raw=False,
        swap_axes=True,
        show=False,
        save=f"_dpt_pseudotime.png")
    
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
    
    if branching>=1:
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
    parser.add_argument('--n_genes', type=int, default = 20)
    parser.add_argument('--clustering_distance', type=float, default = 0.5)
    parser.add_argument('--celltype_key', type=str, required=True)
    parser.add_argument('--root_ctype', type=str, required=True)
    parser.add_argument('--plotdir', type=str, required=True)
    parser.add_argument('--celltype_mask', type=str,  nargs='+', required=True)
    parser.add_argument('--cat_order', type=str,  nargs='+', required=True)
    parser.add_argument('--branching', type=float, default = 0)
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
    n_genes = args.n_genes
    clustering_distance = args.clustering_distance
    cat_order = args.cat_order
    branching = args.branching

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
    adata = compute_pseudotime(adata,celltype_key,root_ctype,diffusion_component,branching)
    stop = timeit.default_timer()
    print(f"Computing diffusion pseudotime in {round(stop-start,2)}s")

    print("3. Identifying dynamic genes")
    start = timeit.default_timer()
    adata_dyn, gene_clusters = identify_dynamic_genes(adata,n_genes,clustering_distance,cat_order,celltype_mask) 
    stop = timeit.default_timer()
    print(f"Dynamic genes identified in {round(stop-start,2)}s")

    print("4. Plotting")
    start = timeit.default_timer()
    plotting(adata,adata_dyn,gene_clusters,celltype_key,branching)
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
