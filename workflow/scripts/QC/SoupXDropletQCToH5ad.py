import argparse
import timeit
import anndata as ad
import scanpy as sc
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import numpy as np
import os


def initialize_parser():
    parser = argparse.ArgumentParser(
        description="Convert SoupX output to h5ad for downstream pipeline"
    )
    parser.add_argument("--raw_h5", type=str, required=True,
                        help="Cell Ranger raw_feature_bc_matrix.h5 (for raw count layer)")
    parser.add_argument("--soupx_dir", type=str, required=True,
                        help="Output directory from SoupX.R")
    parser.add_argument("--output", type=str, required=True, help="Output h5ad path")
    parser.add_argument("--sample", type=str, required=True)
    return parser


def main():
    parser = initialize_parser()
    args = parser.parse_args()

    corrected_dir = os.path.join(args.soupx_dir, "corrected")
    cell_qc_path = os.path.join(args.soupx_dir, "cell_qc.tsv")

    # ------------------------------------------------------------------
    # 1. Load SoupX-corrected counts (10x sparse MTX format)
    # ------------------------------------------------------------------
    print("1. Loading SoupX-corrected matrix")
    start = timeit.default_timer()
    adata = sc.read_10x_mtx(corrected_dir, var_names="gene_symbols", gex_only=True)
    print(f"   Corrected matrix: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"   Loaded in {round(timeit.default_timer()-start, 2)}s")

    # ------------------------------------------------------------------
    # 2. Load raw counts and add as layer
    # ------------------------------------------------------------------
    print("2. Loading raw matrix from Cell Ranger h5")
    start = timeit.default_timer()
    adata_raw = sc.read_10x_h5(args.raw_h5)
    # Handle Multiome h5 (may have ATAC + GEX)
    if "feature_types" in adata_raw.var.columns:
        gex_mask = adata_raw.var["feature_types"] == "Gene Expression"
        adata_raw = adata_raw[:, gex_mask]
    # Align to corrected barcodes and genes
    shared_obs = adata.obs_names.intersection(adata_raw.obs_names)
    shared_var = adata.var_names.intersection(adata_raw.var_names)
    adata_raw_aligned = adata_raw[shared_obs, shared_var]
    # Store corrected in X, raw in layers['raw']
    adata = adata[shared_obs, shared_var].copy()
    adata.layers["raw"] = adata_raw_aligned.X.copy()
    print(f"   Raw counts added as layers['raw'] in {round(timeit.default_timer()-start, 2)}s")

    # ------------------------------------------------------------------
    # 3. Add per-cell QC from SoupX / DropletQC
    # ------------------------------------------------------------------
    print("3. Adding per-cell QC metadata")
    start = timeit.default_timer()
    cell_qc = pd.read_csv(cell_qc_path, sep="\t", index_col="barcode")
    # Align to adata barcodes; warn about any missing barcodes
    missing = adata.obs_names.difference(cell_qc.index)
    if len(missing) > 0:
        print(f"   WARNING: {len(missing)} barcodes in adata not found in cell_qc.tsv; "
              "their QC columns will be NaN.")
    cell_qc = cell_qc.reindex(adata.obs_names)
    for col in cell_qc.columns:
        adata.obs[col] = cell_qc[col].values
    # Add a generic contamination_fraction column for downstream compatibility
    if "soupx_rho" in adata.obs.columns:
        adata.obs["contamination_fraction"] = adata.obs["soupx_rho"]
    print(f"   QC metadata added in {round(timeit.default_timer()-start, 2)}s")

    # ------------------------------------------------------------------
    # 4. Add sample metadata and write
    # ------------------------------------------------------------------
    print("4. Finalising and writing h5ad")
    start = timeit.default_timer()
    adata.obs["sample"] = args.sample
    adata.var_names_make_unique()
    adata.write_h5ad(args.output)
    print(f"   Written to {args.output} in {round(timeit.default_timer()-start, 2)}s")


if __name__ == "__main__":
    start_total = timeit.default_timer()
    main()
    stop_total = timeit.default_timer()
    print(f"Total runtime: {round(stop_total-start_total, 2)}s")
