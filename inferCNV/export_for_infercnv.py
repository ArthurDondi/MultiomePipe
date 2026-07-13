#!/usr/bin/env python
"""
Export a MultiomePipe batch-corrected .h5ad into the flat inputs that inferCNV
(run in R / RStudio) expects.

Run this once on the cluster where the .h5ad lives (in the `scverse` conda env,
which already has anndata / scipy / numpy / pandas), then copy the small output
folder to wherever you run RStudio and point `run_infercnv.R` at it.

inferCNV needs three things (https://github.com/broadinstitute/inferCNV/wiki):
  1. a RAW-counts matrix, genes x cells               -> counts.mtx.gz (+ genes.tsv, barcodes.tsv)
  2. an annotations file, cell -> group               -> annotations_celltype.tsv (and a richer metadata.tsv)
  3. a gene/chromosome ordering file, gene -> chr,start,end -> gene_ordering.tsv (built from the GTF)

We deliberately export the RAW counts from `layers['counts']` (NOT the
log-normalised `X`): inferCNV models raw expression itself.

The matrix is written as a gzipped MatrixMarket file (sparse) so the full
28-sample cohort stays small; `run_infercnv.R` reads it back with
`Matrix::readMM()`.

Example (defaults are filled from config/config_BMO_combined.yaml, so on the
cluster you can usually just run it with no arguments):

    conda activate scverse
    python inferCNV/export_for_infercnv.py \
        --input   /nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/QC/RNA/Merged/BatchCorrection/merged.batch_corrected.h5ad \
        --outdir  /nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/inferCNV/input \
        --gtf     /nobackup/lab_taschner-mandl/arthurdondi/resources/references/hg38/BMO/gencode.v50.basic.annotation.plusGFP.gtf
"""

import os
import gzip
import argparse

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io as sio
import anndata as ad


# Defaults taken from config/config_BMO_combined.yaml so the script is runnable
# with no arguments on the cluster. Override any of them on the command line.
DEF_OUT = "/nobackup/lab_taschner-mandl/arthurdondi/projects/BMO"
DEF_INPUT = f"{DEF_OUT}/QC/RNA/Merged/BatchCorrection/merged.batch_corrected.h5ad"
DEF_OUTDIR = f"{DEF_OUT}/inferCNV/input"
DEF_GTF = ("/nobackup/lab_taschner-mandl/arthurdondi/resources/references/"
           "hg38/BMO/gencode.v50.basic.annotation.plusGFP.gtf")

# Cell-type annotation column, in order of preference. At the BatchCorrection
# stage there is no manual `cell_type` yet, so the scANVI prediction `C_scANVI`
# is normally the best label available; fall back to leiden clusters otherwise.
ANNOTATION_CANDIDATES = ["cell_type", "C_scANVI", "scanvi_labels",
                         "leiden_res_1.00", "leiden_res_0.50"]

# Per-cell metadata columns to carry over (only those present are written).
META_CANDIDATES = ["sample", "donor", "dataset", "barcodes",
                   "leiden_res_0.20", "leiden_res_0.50", "leiden_res_1.00",
                   "C_scANVI", "scanvi_labels", "phase",
                   "total_counts", "n_genes_by_counts", "pct_counts_mt"]

# Chromosomes to keep in the gene-ordering file (drop scaffolds, chrM and the
# GFP transgene contig: none of those carry meaningful chromosomal CNV signal).
STD_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CHROM_RANK = {c: i for i, c in enumerate(STD_CHROMS)}


def pick_annotation_key(adata, requested):
    if requested:
        if requested not in adata.obs.columns:
            raise SystemExit(
                f"[export] requested --annotation-key '{requested}' not in obs. "
                f"Available: {list(adata.obs.columns)}")
        return requested
    for key in ANNOTATION_CANDIDATES:
        if key in adata.obs.columns:
            return key
    raise SystemExit(
        "[export] none of the expected annotation columns "
        f"{ANNOTATION_CANDIDATES} are in obs; pass --annotation-key explicitly. "
        f"Available: {list(adata.obs.columns)}")


def get_raw_counts(adata, layer):
    """genes x cells sparse matrix of RAW counts (transpose of AnnData's X)."""
    if layer in adata.layers:
        X = adata.layers[layer]
    else:
        print(f"[export] WARNING: layer '{layer}' absent; falling back to .X. "
              "If .X is log-normalised this is NOT what inferCNV wants — "
              "re-export once a raw-counts layer is available.")
        X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    # genes (rows) x cells (cols)
    mat = X.T.tocsc()
    # Counts should be integral; store as integers so the .mtx is compact and
    # unambiguous. If they are not (unexpected), keep floats and warn.
    data = mat.data
    if data.size and np.allclose(data, np.rint(data)):
        mat.data = np.rint(data).astype(np.int32)
    else:
        print("[export] WARNING: counts layer is not integer-valued; writing as float.")
    return mat


def write_matrix(mat, genes, cells, outdir):
    mtx_path = os.path.join(outdir, "counts.mtx.gz")
    print(f"[export] writing {mat.shape[0]} genes x {mat.shape[1]} cells -> {mtx_path}")
    with gzip.open(mtx_path, "wb") as fh:
        sio.mmwrite(fh, mat, comment="genes x cells raw counts (inferCNV input)")
    # One entry per line, no header: read back in R with readLines().
    pd.Series(genes).to_csv(os.path.join(outdir, "genes.tsv"),
                            index=False, header=False)
    pd.Series(cells).to_csv(os.path.join(outdir, "barcodes.tsv"),
                            index=False, header=False)


def write_annotations(adata, ann_key, outdir):
    cells = adata.obs_names
    labels = adata.obs[ann_key].astype(str).fillna("Unknown")

    # 2-column, no header: cell <TAB> group. This is the direct inferCNV
    # `annotations_file` for a joint (all-samples) run.
    ann = pd.DataFrame({"cell": cells, "group": labels.values})
    ann.to_csv(os.path.join(outdir, "annotations_celltype.tsv"),
               sep="\t", index=False, header=False)

    # Richer per-cell metadata (with header) so run_infercnv.R can build
    # per-sample annotations and pick reference groups.
    meta = pd.DataFrame(index=cells)
    meta.insert(0, "cell", cells)
    meta["celltype"] = labels.values
    for col in META_CANDIDATES:
        if col in adata.obs.columns and col not in meta.columns:
            meta[col] = adata.obs[col].astype(str).fillna("NA").values
    meta.to_csv(os.path.join(outdir, "metadata.tsv"), sep="\t", index=False)

    print(f"[export] annotation key = '{ann_key}'")
    print("[export] cells per group:")
    print(labels.value_counts().to_string())
    if "sample" in adata.obs.columns:
        print(f"[export] samples = {adata.obs['sample'].nunique()}")
    return meta


def _open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def write_gene_ordering(gtf_path, genes, outdir):
    """Build gene_ordering.tsv (gene <TAB> chr <TAB> start <TAB> end) from a GTF.

    Only `gene` features on standard chromosomes are kept; the file is limited
    to genes present in the matrix and sorted in genomic order (inferCNV walks
    genes in file order along each chromosome).
    """
    if not gtf_path or not os.path.exists(gtf_path):
        print(f"[export] GTF not found ({gtf_path}); SKIPPING gene_ordering.tsv. "
              "Provide --gtf, or supply inferCNV's stock gencode positions file "
              "to run_infercnv.R.")
        return

    import re
    name_re = re.compile(r'gene_name "([^"]+)"')
    rows = {}  # gene_name -> [chrom, start, end]
    with _open_maybe_gzip(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom = f[0]
            if chrom not in CHROM_RANK:
                continue
            m = name_re.search(f[8])
            if not m:
                continue
            name = m.group(1)
            start, end = int(f[3]), int(f[4])
            if name in rows:  # same symbol seen twice -> widen the span
                prev = rows[name]
                if prev[0] == chrom:
                    prev[1] = min(prev[1], start)
                    prev[2] = max(prev[2], end)
            else:
                rows[name] = [chrom, start, end]

    gene_set = set(map(str, genes))
    kept = [(g, c, s, e) for g, (c, s, e) in rows.items() if g in gene_set]
    kept.sort(key=lambda r: (CHROM_RANK[r[1]], r[2], r[3]))

    out = os.path.join(outdir, "gene_ordering.tsv")
    with open(out, "w") as fh:
        for g, c, s, e in kept:
            fh.write(f"{g}\t{c}\t{s}\t{e}\n")
    n_matrix = len(gene_set)
    print(f"[export] gene ordering: {len(kept)}/{n_matrix} matrix genes "
          f"placed on chr1-22/X/Y -> {out}")
    missing = n_matrix - len(kept)
    if missing:
        print(f"[export] {missing} matrix genes have no standard-chromosome "
              "position (scaffolds / GFP / unmatched symbols); inferCNV ignores them.")


def main():
    p = argparse.ArgumentParser(
        description="Export a batch-corrected .h5ad into inferCNV inputs.")
    p.add_argument("--input", default=DEF_INPUT, help="batch-corrected .h5ad")
    p.add_argument("--outdir", default=DEF_OUTDIR, help="output directory")
    p.add_argument("--gtf", default=DEF_GTF,
                   help="GTF used to build gene_ordering.tsv (pass '' to skip)")
    p.add_argument("--counts-layer", default="counts",
                   help="AnnData layer holding raw counts (default: counts)")
    p.add_argument("--annotation-key", default=None,
                   help=f"obs column for cell groups (default: first of "
                        f"{ANNOTATION_CANDIDATES} present)")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"[export] reading {args.input}")
    adata = ad.read_h5ad(args.input)
    print(f"[export] {adata.n_obs} cells x {adata.n_vars} genes")

    ann_key = pick_annotation_key(adata, args.annotation_key)

    mat = get_raw_counts(adata, args.counts_layer)
    write_matrix(mat, list(adata.var_names), list(adata.obs_names), args.outdir)
    write_annotations(adata, ann_key, args.outdir)
    write_gene_ordering(args.gtf, adata.var_names, args.outdir)

    print(f"[export] DONE -> {args.outdir}")
    print("[export] next: run inferCNV/run_infercnv.R in RStudio, pointing "
          "INPUT_DIR at this folder.")


if __name__ == "__main__":
    main()
