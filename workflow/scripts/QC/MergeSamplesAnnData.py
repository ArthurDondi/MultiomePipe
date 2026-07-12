import os
import argparse
import timeit
import colorsys
import anndata as ad
import scanpy as sc
import pandas as pd
import json

# -----------------------------
# Categorical palette (large N)
# -----------------------------
def build_distinct_palette(n):
    """Return `n` visually distinct hex colours for a categorical variable.

    scanpy falls back to low-contrast / shaded palettes beyond ~28 categories,
    and because samples are ordered by name, look-alike colours end up next to
    each other. We walk the HSV wheel in golden-angle (~137.5 deg) steps so
    *consecutive* categories (alphabetically adjacent sample names) land far
    apart in hue, and rotate through four saturation/value tiers so any near-hue
    repeat at high `n` still differs in brightness. Pure stdlib — no new deps.
    """
    golden = 0.6180339887498949  # 1/phi -> ~137.5 deg hue steps
    tiers = ((0.78, 1.00), (1.00, 0.62), (0.48, 0.95), (0.95, 0.42))
    palette = []
    for i in range(n):
        h = (i * golden) % 1.0
        s, v = tiers[i % len(tiers)]
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        palette.append("#{:02x}{:02x}{:02x}".format(
            round(r * 255), round(g * 255), round(b * 255)))
    return palette


def apply_distinct_palette(adata, keys, min_categories=21):
    """Attach the high-contrast palette to big categorical obs columns.

    Only columns with at least `min_categories` levels are overridden — smaller
    ones keep scanpy's already-distinct tab10 / default_20. Colours are written
    to ``uns['<key>_colors']``, which scanpy honours for every subsequent plot
    and which persists into the written h5ad.
    """
    for key in keys:
        if not key or key not in adata.obs.columns:
            continue
        col = adata.obs[key]
        if not isinstance(col.dtype, pd.CategoricalDtype):
            col = col.astype("category")
            adata.obs[key] = col
        n = len(col.cat.categories)
        if n >= min_categories:
            adata.uns[f"{key}_colors"] = build_distinct_palette(n)


def _sanitize_for_write(adata):
    """Make an outer-joined AnnData writable to .h5ad.

    Two artefacts of the union (outer) merge trip anndata's h5ad writer:

    * The var index is *named* 'gene_symbols' while a differing 'gene_symbols'
      column exists (var_names_make_unique suffixes duplicate symbols such as
      'TBCE' -> 'TBCE-1'); anndata refuses to write an index whose name matches
      a column with different values.
    * Boolean columns (the mt / ribo / hb QC flags) get upcast to ``object``
      dtype when reindexing to the gene union introduces NaN. anndata then
      serialises them as strings and h5py raises "Can't implicitly convert
      non-string objects to strings".

    Clear any index name that collides with a column, and cast all-boolean
    object columns back to a real bool dtype (missing -> False). String and
    numeric columns are left untouched.
    """
    for frame in (adata.var, adata.obs):
        if frame.index.name in frame.columns:
            frame.index.name = None
        for col in frame.columns:
            s = frame[col]
            if s.dtype != object:
                continue
            non_null = s.dropna()
            if len(non_null) and all(
                type(v).__name__ in ("bool", "bool_") for v in non_null
            ):
                frame[col] = s.fillna(False).astype(bool)


# -----------------------------
# Functions
# -----------------------------
def merge_data(input_files, samples, donor_key, sample_key, dataset_key=None):

    # `input_files[i]` and `samples[i]` arrive parallel: the MergeSamplesAnnData
    # rule builds both `--input` and `--samples` from the same SAMPLES order, so
    # file i already belongs to sample i. Sorting the two lists *independently*
    # (the previous `sorted(input_files)` / `sorted(samples)`) silently re-pairs a
    # file with the WRONG sample label whenever one sample name is a lexicographic
    # prefix of another followed by a character that sorts before '/'. Each file is
    # `QC/RNA/<sample>/Filtering/<sample>_filtered.h5ad`, so after the sample dir
    # comes '/' (0x2f); '-' (0x2d) and '.' (0x2e) sort *before* it. So e.g.
    # samples 'BMO' + 'BMO-IMR-DOX' (this project's own hyphen scheme) or
    # 'S1' + 'S1.rerun' sort their file paths in a different order than their
    # names, and every affected sample's cells get another sample's
    # name/donor/dataset — a silent label swap that corrupts all downstream batch
    # covariates and per-sample biology. Sort the (file, sample) PAIRS together so
    # the ordering stays deterministic (sorted by sample, as before) *and* the
    # file->label pairing can never diverge.
    paired = sorted(zip(input_files, samples), key=lambda fs: fs[1])
    samples = [sample for _, sample in paired]

    adatas = []
    for f, sample in paired:
        adata = ad.io.read_h5ad(f)
        adata.var_names_make_unique()
        adata.obs['barcodes'] = adata.obs_names
        adata.obs_names = adata.obs_names + '-' + sample
        adatas.append(adata)
    print("Concat")
    # Outer join keeps the UNION of genes across samples (missing entries filled
    # with 0), rather than the intersection. An inner join silently dropped any
    # gene not detected in *every* sample — most notably transgenes such as EGFP,
    # which are only expressed in a subset of samples and so vanished from the
    # merged object, breaking downstream `adata[:, "EGFP"]` / colour lookups.
    adata_concat = ad.concat(
        adatas,
        merge="first",
        uns_merge="first",
        join="outer",
        fill_value=0,
        label="sample",
        # `samples` is now ordered to match `adatas` (both come from `paired`),
        # so block i is correctly labelled with sample i.
        keys=samples
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
    # High-contrast palette for the many-category batch covariates (samples /
    # donors), so alphabetically adjacent samples are not drawn in look-alike
    # shades. Small columns (e.g. dataset) keep scanpy's default palette.
    apply_distinct_palette(adata, [dataset_key, donor_key, sample_key])
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
    _sanitize_for_write(adata)
    adata.write_h5ad(output_file)
    stop = timeit.default_timer()
    print(f"Data written in {round(stop-start,2)}s")

if __name__ == "__main__":
    start_total = timeit.default_timer()
    main()
    stop_total = timeit.default_timer()
    print(f"Total runtime: {round(stop_total-start_total,2)}s")
