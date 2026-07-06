import os
import argparse
import timeit
import anndata as ad
import scanpy as sc
import pandas as pd
import json

# -----------------------------
# Functions
# -----------------------------
def merge_data(input_files, samples, donor_key, sample_key, dataset_key=None):

    adatas = [ad.io.read_h5ad(f) for f in sorted(input_files)]
    #Adding sample name to cell names — use sorted(samples) to match sorted file order
    for adata, sample in zip(adatas, sorted(samples)):
        adata.var_names_make_unique()
        adata.obs['barcodes'] = adata.obs_names
        adata.obs_names = adata.obs_names+'-'+sample
    print("Concat")
    adata_concat = ad.concat(
        adatas,
        merge="first",
        uns_merge="first",
        join="inner",
        label="sample",
        keys=sorted(samples)
    )
    print("Concat done")
    del adatas
    # The cellbender background-correction path stores corrected counts in a
    # 'cellbender' layer; the soupx path already has corrected counts in X.
    if 'cellbender' in adata_concat.layers:
        adata_concat.X = adata_concat.layers['cellbender']

    # Ensure batch columns use a consistent sorted categorical order
    # so that colour palettes are assigned identically in every UMAP panel.
    for col in [donor_key, sample_key, dataset_key]:
        if col is None:
            continue
        if col in adata_concat.obs.columns:
            adata_concat.obs[col] = pd.Categorical(
                adata_concat.obs[col],
                categories=sorted(adata_concat.obs[col].unique())
            )

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


def run_normalization_and_clustering(adata, donor_key, sample_key, dataset_key, is_filtered):

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000,batch_key=sample_key)
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
    adata.obsm["X_umap_nocorrection"] = adata.obsm["X_umap"].copy()
    for res in [0.2, 0.5, 1]:
        sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res) # cluster labels in adata.obs['leiden']
    sc.pl.umap(adata,
            color=[f"leiden_res_{res:4.2f}" for res in [0.2, 0.5, 1]],
            legend_loc="on data",
            show=False,
            save=f"_leiden_res_merge.png")
    
    QCs = ["log1p_total_counts", "pct_counts_mt", "doublet_score"]
    if not is_filtered:
        # Nuclear fraction (DropletQC, SoupX path): colour the UMAP by nf so the
        # high-nf damaged population is visible amongst the other QC metrics,
        # mirroring the per-sample QC UMAP in RawFilteringRNA.
        if "nf_umi" in adata.obs.columns and adata.obs["nf_umi"].notna().any():
            QCs.append("nf_umi")
        # The ambient-contamination metric is named differently depending on the
        # background-correction backend: CellBender writes 'background_fraction',
        # SoupX writes 'soupx_rho' (aliased to 'contamination_fraction'). Overlay
        # whichever one is present so the UMAP does not crash with a KeyError.
        for col in ["background_fraction", "contamination_fraction", "soupx_rho"]:
            if col in adata.obs.columns:
                QCs.append(col)
                break
    # Guard against any QC column that did not survive the merge/join.
    QCs = [c for c in QCs if c in adata.obs.columns]
    sc.pl.umap(adata,
                color=QCs,
                show=False,
                save=f"_QC_merge.png")
    # Colour the uncorrected UMAP by dataset-of-origin (the coarsest batch
    # covariate) and sample. `wspace` widens the gap between the two panels so
    # the left panel's legend is not overpainted by the right panel; scanpy
    # saves with bbox_inches='tight', so the right panel's (long) sample legend
    # is fully included. Fall back to donor if no dataset column is present.
    batch_key = dataset_key if (dataset_key and dataset_key in adata.obs.columns) else donor_key
    sc.pl.umap(adata,
                color=[batch_key,sample_key],
                wspace=0.5,
                show=False,
                save=f"_batch_nocorrection_merge.png")
    
    return adata

def initialize_parser():
    parser = argparse.ArgumentParser(description='Merging adatas from all samples and plotting')
    parser.add_argument('--input', type=str, nargs='+', required=True, help="space separated h5ad list")
    parser.add_argument('--output', type=str, required=True, help="h5ad")
    parser.add_argument('--samples', type=str, nargs='+', required=True, help="space separated sample names list")
    parser.add_argument('--plotdir', type=str, required=True)
    parser.add_argument('--sample_key', type=str, required=True)
    parser.add_argument('--donor_key', type=str, required=True)
    parser.add_argument('--dataset_key', type=str, default=None)
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
    samples = args.samples
    plotdir = args.plotdir
    donor_key = args.donor_key
    sample_key = args.sample_key
    dataset_key = args.dataset_key
    is_filtered = args.is_filtered

    sc.settings.figdir = plotdir

    # 1. Merge data
    print("1. Merge data")
    start = timeit.default_timer()
    adata = merge_data(input_files, samples, donor_key, sample_key, dataset_key)
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
    adata = run_normalization_and_clustering(adata, donor_key, sample_key, dataset_key, is_filtered)
    stop = timeit.default_timer()
    print(f"Normalization + clustering done in {round(stop-start,2)}s")

    # 4. Write data
    print("4. Write data")
    start = timeit.default_timer()
    adata.write_h5ad(output_file)
    stop = timeit.default_timer()
    print(f"Data written in {round(stop-start,2)}s")

if __name__ == "__main__":
    start_total = timeit.default_timer()
    main()
    stop_total = timeit.default_timer()
    print(f"Total runtime: {round(stop_total-start_total,2)}s")
