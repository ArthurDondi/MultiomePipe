import argparse
import timeit
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd


def get_raw_counts(adata):
    """Get raw counts from adata, checking layers in order of priority."""
    for layer in ['counts', 'cellbender']:
        if layer in adata.layers:
            return adata.layers[layer].copy()
    return adata.X.copy()


def align_genes(adata_query, adata_ref):
    """Align genes between query and reference adata.

    Tries direct var_names intersection first, then falls back to matching
    via 'feature_name' or 'gene_names' columns in the reference var table.

    Returns (adata_query_sub, adata_ref_sub) subsetted to common genes with
    aligned var_names.
    """
    common = list(adata_query.var_names.intersection(adata_ref.var_names))

    if len(common) >= 100:
        return adata_query[:, common].copy(), adata_ref[:, common].copy()

    # Fall back to matching query gene symbols against reference gene-name columns
    for ref_col in ['feature_name', 'gene_names']:
        if ref_col not in adata_ref.var.columns:
            continue
        gene_map = dict(zip(adata_ref.var[ref_col], adata_ref.var_names))
        matched_query = [g for g in adata_query.var_names if g in gene_map]
        if len(matched_query) < 100:
            continue
        matched_ref = [gene_map[g] for g in matched_query]
        adata_ref_sub = adata_ref[:, matched_ref].copy()
        adata_ref_sub.var_names = matched_query
        adata_query_sub = adata_query[:, matched_query].copy()
        print(f"Gene alignment via reference column '{ref_col}': {len(matched_query)} genes")
        return adata_query_sub, adata_ref_sub

    raise ValueError(
        f"Unable to align genes: only {len(common)} direct var_names matches found "
        "and no usable gene-name column ('feature_name', 'gene_names') in reference."
    )


def label_transfer(adata, reference_h5ad, columns, celltype_key, n_neighbors, plotdir):
    """Transfer labels from a reference dataset to the query using scanpy ingest.

    Parameters
    ----------
    adata : AnnData
        Query dataset (batch-corrected merged h5ad).
    reference_h5ad : str
        Path to the reference h5ad file.
    columns : list of str
        Columns in reference adata.obs whose labels should be transferred.
    celltype_key : str
        Primary cell-type key to set in adata.obs (set to columns[0]).
    n_neighbors : int
        Number of neighbours for KNN classification.
    plotdir : str
        Directory for output plots.
    """
    print("  Loading reference dataset ...")
    adata_ref = ad.io.read_h5ad(reference_h5ad)

    # Validate requested columns
    missing_cols = [c for c in columns if c not in adata_ref.obs.columns]
    if missing_cols:
        raise ValueError(f"Columns not found in reference obs: {missing_cols}")

    # Align genes
    adata_q_sub, adata_r_sub = align_genes(adata, adata_ref)
    n_common = adata_q_sub.n_vars
    print(f"  Common genes: {n_common}")

    # ---- Process reference -----------------------------------------------
    # Use raw counts when available, otherwise use X as-is
    adata_r_sub.X = get_raw_counts(adata_r_sub)
    sc.pp.normalize_total(adata_r_sub, target_sum=1e4)
    sc.pp.log1p(adata_r_sub)
    n_hvg = min(2000, n_common - 1)
    sc.pp.highly_variable_genes(adata_r_sub, n_top_genes=n_hvg)
    hvg_mask = adata_r_sub.var['highly_variable']
    hvg_genes = adata_r_sub.var_names[hvg_mask].tolist()

    adata_r_hvg = adata_r_sub[:, hvg_genes].copy()
    sc.pp.scale(adata_r_hvg, max_value=10)
    sc.pp.pca(adata_r_hvg)
    sc.pp.neighbors(adata_r_hvg, n_neighbors=n_neighbors)

    # ---- Process query on the same HVGs ----------------------------------
    query_hvg_genes = [g for g in hvg_genes if g in adata_q_sub.var_names]
    if len(query_hvg_genes) < 50:
        raise ValueError(
            f"Too few HVG genes in query ({len(query_hvg_genes)}). "
            "Cannot perform label transfer."
        )

    adata_q_hvg = adata_q_sub[:, query_hvg_genes].copy()
    adata_q_hvg.X = get_raw_counts(adata_q_hvg)
    sc.pp.normalize_total(adata_q_hvg, target_sum=1e4)
    sc.pp.log1p(adata_q_hvg)
    sc.pp.scale(adata_q_hvg, max_value=10)

    # If any HVG genes were missing in the query, re-subset reference to match
    if len(query_hvg_genes) < len(hvg_genes):
        print(
            f"  {len(hvg_genes) - len(query_hvg_genes)} HVG genes absent in query; "
            "re-subsetting reference."
        )
        adata_r_hvg = adata_r_hvg[:, query_hvg_genes].copy()
        sc.pp.pca(adata_r_hvg)
        sc.pp.neighbors(adata_r_hvg, n_neighbors=n_neighbors)

    # ---- Ingest (KNN label transfer) -------------------------------------
    sc.tl.ingest(adata_q_hvg, adata_r_hvg, obs=columns)

    # ---- Copy transferred labels back to query adata ---------------------
    for col in columns:
        adata.obs[col] = adata_q_hvg.obs[col].values
        print(f"  Transferred '{col}': {adata.obs[col].nunique()} unique labels")

    # Set primary celltype_key to the first transferred column
    adata.obs[celltype_key] = adata.obs[columns[0]].values

    # ---- Plots -----------------------------------------------------------
    sc.settings.figdir = plotdir
    for col in columns:
        sc.pl.umap(
            adata,
            color=col,
            legend_loc="on data",
            legend_fontsize=8,
            show=False,
            save=f"_label_transfer_{col}.png"
        )

    return adata


def initialize_parser():
    parser = argparse.ArgumentParser(
        description='Label transfer from reference h5ad to query h5ad using scanpy ingest'
    )
    parser.add_argument('--input', type=str, required=True, help="Query h5ad (batch-corrected)")
    parser.add_argument('--output', type=str, required=True, help="Output h5ad with transferred labels")
    parser.add_argument('--reference', type=str, required=True, help="Reference h5ad")
    parser.add_argument(
        '--reference_columns', type=str, nargs='+', required=True,
        help="Reference obs columns to transfer (space-separated)"
    )
    parser.add_argument('--celltype_key', type=str, default='cell_type',
                        help="Primary cell-type key to set in obs (default: 'cell_type')")
    parser.add_argument('--n_neighbors', type=int, default=15,
                        help="Number of neighbours for KNN (default: 15)")
    parser.add_argument('--plotdir', type=str, required=True)
    return parser


# -----------------------------
# Main
# -----------------------------
def main():
    parser = initialize_parser()
    args = parser.parse_args()

    sc.settings.figdir = args.plotdir

    print("1. Load query data")
    start = timeit.default_timer()
    adata = ad.io.read_h5ad(args.input)
    print(f"   Query: {adata}")
    stop = timeit.default_timer()
    print(f"   Loaded in {round(stop - start, 2)}s")

    print("2. Label transfer")
    start = timeit.default_timer()
    adata = label_transfer(
        adata,
        reference_h5ad=args.reference,
        columns=args.reference_columns,
        celltype_key=args.celltype_key,
        n_neighbors=args.n_neighbors,
        plotdir=args.plotdir
    )
    stop = timeit.default_timer()
    print(f"   Label transfer done in {round(stop - start, 2)}s")

    print("3. Write data")
    start = timeit.default_timer()
    adata.write_h5ad(args.output)
    stop = timeit.default_timer()
    print(f"   Written in {round(stop - start, 2)}s")


if __name__ == "__main__":
    start_total = timeit.default_timer()
    main()
    stop_total = timeit.default_timer()
    print(f"Total runtime: {round(stop_total - start_total, 2)}s")
