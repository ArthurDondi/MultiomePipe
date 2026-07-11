import os
import re
import json
import argparse
import timeit
import colorsys
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import pandas as pd


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


def read_h5ad(path):
    """Read an h5ad, tolerating ``null``-encoded scalars from a newer anndata.

    ``sc.pp.log1p`` records ``uns['log1p'] = {'base': None}``. anndata >= 0.12
    (the version in envs/scverse.yaml, which writes the merged object) serialises
    that ``None`` with ``encoding_type='null'``. The 0.11.x anndata pinned in
    envs/scvi.yaml — used by this script for the scvi/scanvi backends — has no
    read method registered for ``null`` and blows up with an ``IORegistryError``.

    A plain ``read_dispatched`` callback can't recover, because anndata resolves
    the (missing) read method *before* the callback runs. So we intercept at the
    enclosing mapping: whenever a dict directly contains a ``null``-encoded child,
    we read that dict ourselves and substitute ``None``, exactly the value the
    newer anndata would have produced. Everything else defers to the default
    reader, so behaviour is unchanged on anndata versions that read ``null`` fine.
    """
    try:
        return ad.io.read_h5ad(path)
    except Exception as err:
        if "encoding_type='null'" not in str(err):
            raise

        import h5py
        from anndata.experimental import read_dispatched

        def _is_null(elem):
            return elem.attrs.get("encoding-type") == "null"

        def callback(read_func, elem_name, elem, iospec):
            if iospec.encoding_type == "null":
                return None
            if iospec.encoding_type == "dict" and any(
                _is_null(v) for v in elem.values()
            ):
                return {
                    k: None if _is_null(v) else read_dispatched(v, callback)
                    for k, v in elem.items()
                }
            return read_func(elem)

        try:
            with h5py.File(path, "r") as f:
                return read_dispatched(f, callback)
        except Exception:
            # Fall back to the original error so genuine read problems aren't
            # masked by the null-encoding workaround.
            raise err


def prepare_annotation(adata, markers_file):

    with open(markers_file, "r") as f:
        marker_genes = json.load(f)

    # Keeping only genes in dataset
    marker_genes_in_data = {}
    for ct, markers in marker_genes.items():
        markers_found = []
        for marker in markers:
            if marker in adata.var.index:
                markers_found.append(marker)
        marker_genes_in_data[ct] = markers_found

    for res in [0.2, 0.5, 1]:
        sc.pl.dotplot(
            adata,
            groupby=f"leiden_res_{res:4.2f}",
            var_names=marker_genes_in_data,
            standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
            show=False,
            save=f"_manual_markers_leiden_res_{res:4.2f}_merge.png"
        )

    return adata


def _plot_umap_on_data_with_legend(adata, col, unlabeled_category=None, save=None):
    """UMAP with small, non-bold on-data labels *and* a right-margin legend.

    scanpy's `legend_loc` is either 'on data' or 'right margin', never both, so
    we draw the small on-data labels through scanpy and add the categorical
    colour legend ourselves from the palette scanpy stores in `uns`. Used for the
    scANVI label / prediction overlays, which have many fine-grained cell types
    whose on-data labels overlap — the side legend stays readable regardless.

    If `unlabeled_category` is given and present in `col`, those seed cells are
    forced to light grey so the reference-matched labels stand out (only the
    scanvi_labels overlay carries it; C_scANVI predictions never do).
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    # Ensure a stable categorical so the palette scanpy writes to uns lines up
    # with the categories we build the legend from (scANVI predictions come back
    # as a plain object column, not a categorical).
    if not isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
        adata.obs[col] = adata.obs[col].astype("category")
    cats = list(adata.obs[col].cat.categories)

    # Force the unlabeled/"Unknown" seed cells to light grey so the reference-
    # matched labels stand out. We pre-set the palette in uns — giving every
    # other category the same default colour scanpy would have picked (same
    # size thresholds) — so both the scatter and the legend below honour it.
    if unlabeled_category is not None and unlabeled_category in cats:
        n = len(cats)
        base = (sc.pl.palettes.default_20 if n <= 20
                else sc.pl.palettes.default_28 if n <= 28
                else sc.pl.palettes.default_102)
        palette = [base[i % len(base)] for i in range(n)]
        palette[cats.index(unlabeled_category)] = "lightgray"
        adata.uns[f"{col}_colors"] = palette

    fig, ax = plt.subplots(figsize=(8, 7))
    sc.pl.umap(adata,
               color=col,
               ax=ax,
               legend_loc="on data",
               legend_fontsize=5,
               legend_fontweight="normal",
               legend_fontoutline=1,
               show=False)

    colors = adata.uns.get(f"{col}_colors")
    if colors is not None and len(colors) >= len(cats):
        handles = [Patch(facecolor=colors[i], label=str(cats[i]))
                   for i in range(len(cats))]
        ax.legend(handles=handles,
                  loc="center left",
                  bbox_to_anchor=(1.02, 0.5),
                  frameon=False,
                  fontsize=6,
                  ncol=2 if len(cats) > 20 else 1,
                  title=col,
                  title_fontsize=7)

    figdir = str(sc.settings.figdir)
    os.makedirs(figdir, exist_ok=True)
    # Match scanpy's own naming (umap<save>) so downstream expectations hold;
    # `save` overrides the filename for the per-stage diagnostic embeddings.
    fig.savefig(os.path.join(figdir, save or f"umap_{col}_merge.png"),
                dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_stage_embedding(adata, rep_key, label_col, tag, unlabeled_category=None):
    """Diagnostic: UMAP on `rep_key` coloured by the *reference* labels
    `label_col` (no scANVI propagation), saved as umap_stage_<tag>_merge.png.

    Lets you compare how each stage of Method A — scVI, scANVI #1 (integration),
    scANVI #2 (annotation) — reshapes the manifold, all coloured by the same
    input labels. The UMAP coords are stashed under obsm['X_umap_<tag>'] so the
    main clustering UMAP (built later on the integration embedding) is unaffected.
    """
    if rep_key not in adata.obsm or label_col not in adata.obs.columns:
        print(f"  [stage plot] skip '{tag}': missing {rep_key} or {label_col}")
        return
    sc.pp.neighbors(adata, use_rep=rep_key)
    sc.tl.umap(adata)
    adata.obsm[f"X_umap_{tag}"] = adata.obsm["X_umap"].copy()
    _plot_umap_on_data_with_legend(adata, label_col, unlabeled_category,
                                   save=f"umap_stage_{tag}_merge.png")


def clustering(adata, sample_key, donor_key, dataset_key, is_filtered, use_rep,
               extra_colors=None, unlabeled_category=None):
    """Build neighbours/UMAP on `use_rep`, then Leiden-cluster and plot.

    `use_rep` is the batch-corrected embedding produced by the chosen method
    (X_pca_harmony / X_scVI / X_scANVI). `dataset_key` colours the batch UMAP by
    dataset-of-origin (falling back to `donor_key` when absent). `extra_colors`
    are optional extra obs columns to overlay on the UMAP (e.g. scANVI labels /
    predictions); `unlabeled_category` greys those seed cells out in the
    scANVI-label overlay.
    """

    if use_rep in adata.obsm:
        sc.pp.neighbors(adata, use_rep=use_rep)
        sc.tl.umap(adata)

    for res in [0.2, 0.5, 1]:
        sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res)
        sc.tl.rank_genes_groups(adata, groupby=f"leiden_res_{res:4.2f}", method="wilcoxon")
        sc.pl.rank_genes_groups_dotplot(adata,
                                        groupby=f"leiden_res_{res:4.2f}",
                                        standard_scale="var",
                                        n_genes=5,
                                        show=False,
                                        save=f"_diff_markers_leiden_res_{res:4.2f}_merge.png")
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

    # Colour the integration UMAP by dataset-of-origin (the coarsest batch
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
                color=[batch_key, sample_key],
                wspace=0.5,
                show=False,
                save=f"_merge.png")

    for col in (extra_colors or []):
        if col in adata.obs.columns:
            _plot_umap_on_data_with_legend(adata, col, unlabeled_category)


# -----------------------------
# Batch correction backends
# -----------------------------
def run_harmony(adata, batch_keys):
    """Harmony integration -> adata.obsm['X_pca_harmony']."""
    if any(len(adata.obs[k].unique()) > 1 for k in batch_keys):
        sce.pp.harmony_integrate(adata,
                                key=batch_keys,
                                theta=3.0,           # default ~2; smaller = weaker integration
                                sigma=0.1)
    return "X_pca_harmony"


def _combined_batch(adata, batch_keys):
    """Single categorical batch column combining all batch covariates.

    Mirrors Harmony's behaviour of using dataset+donor+sample jointly.
    """
    combined = adata.obs[batch_keys].astype(str).agg("|".join, axis=1)
    adata.obs["_scvi_batch"] = pd.Categorical(combined)
    return "_scvi_batch"


def _select_hvg(adata, n_hvg, batch_key):
    """seurat_v3 HVG on raw counts (layer 'counts'); returns gene subset adata."""
    layer = "counts" if "counts" in adata.layers else None
    n_hvg = min(n_hvg, adata.n_vars)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_hvg,
        flavor="seurat_v3",
        layer=layer,
        batch_key=batch_key,
        subset=False,
    )
    return adata[:, adata.var["highly_variable"]].copy()


def run_scvi(adata, batch_keys, n_hvg, n_latent, max_epochs):
    """Train scVI on raw counts -> adata.obsm['X_scVI']. Returns the model."""
    import scvi
    scvi.settings.seed = 0

    batch_key = _combined_batch(adata, batch_keys)
    adata_hvg = _select_hvg(adata, n_hvg, batch_key)
    if "counts" not in adata_hvg.layers:
        # RawFilteringRNA guarantees a 'counts' layer; fall back to X otherwise.
        adata_hvg.layers["counts"] = adata_hvg.X.copy()

    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        layer="counts",
        batch_key=batch_key,
    )
    model = scvi.model.SCVI(adata_hvg, n_latent=n_latent)
    model.train(max_epochs=max_epochs, accelerator="auto")

    adata.obsm["X_scVI"] = model.get_latent_representation()
    # Keep the HVG-subset handle for scANVI (shares the same setup).
    model.adata_hvg = adata_hvg
    return model


def run_scanvi(adata, scvi_model, labels_key, unlabeled_category, max_epochs,
               rep_key="X_scANVI", pred_key="C_scANVI"):
    """Fine-tune scANVI from a trained scVI model.

    Stores the latent representation in ``adata.obsm[rep_key]`` (when not None)
    and the predicted labels in ``adata.obs[pred_key]`` (when not None). Method A
    calls this twice on one scVI base: once with the coarse *integration* labels
    (keep the embedding used downstream) and once with the *fine* labels (keep the
    predictions used for annotation).
    """
    import scvi

    adata_hvg = scvi_model.adata_hvg
    adata_hvg.obs[labels_key] = adata.obs[labels_key].values

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        adata=adata_hvg,
        labels_key=labels_key,
        unlabeled_category=unlabeled_category,
    )
    scanvi_model.train(max_epochs=max_epochs, accelerator="auto")

    if rep_key is not None:
        adata.obsm[rep_key] = scanvi_model.get_latent_representation()
    if pred_key is not None:
        adata.obs[pred_key] = scanvi_model.predict(adata_hvg)
    return scanvi_model


# -----------------------------
# Label ingestion (for scANVI)
# -----------------------------
def _norm_barcode(bc):
    """Extract the 16bp 10x cell barcode from an arbitrary barcode string.

    Reference h5ads often prefix/suffix barcodes (e.g. 'BMO1_AAAC...-1' or
    'AAAC...-1-2'); we match on the bare nucleotide run so naming differences
    between the published object and our re-quantification don't matter.
    """
    m = re.search(r"[ACGTN]{14,}", str(bc).upper())
    return m.group(0) if m else None


def _apply_ref_column(ref, col, ref_bc, ds_bc, target, dataset, unlabeled_category):
    """Barcode-match one reference column onto `target` for one dataset's cells."""
    if col is None:
        return
    if col not in ref.obs.columns:
        print(f"  [labels] dataset '{dataset}': column '{col}' absent in reference "
              f"obs ({list(ref.obs.columns)[:10]} ...); leaving those cells "
              f"as '{unlabeled_category}'")
        return
    ref_map = pd.Series(ref.obs[col].astype(str).values, index=ref_bc.values)
    ref_map = ref_map[~ref_map.index.duplicated(keep="first")]
    mapped = ds_bc.map(ref_map)
    matched = mapped.notna()
    target.loc[matched[matched].index] = mapped[matched].values
    print(f"  [labels] dataset '{dataset}' [{col}]: matched "
          f"{int(matched.sum())}/{len(ds_bc)} cells")


def build_scanvi_labels(adata, dataset_key, label_datasets, label_refs,
                        label_columns, label_columns_fine, unlabeled_category,
                        int_key="scanvi_labels", fine_key="scanvi_labels_fine"):
    """Assemble TWO label columns from per-dataset reference h5ads (Method A).

    `label_columns` gives the coarse *integration* label per source (drives
    scANVI #1 / the embedding); `label_columns_fine` gives the *fine* annotation
    label (drives scANVI #2 / the predictions). Both are read in a single pass
    over each reference and matched by cell barcode. Cells with no match — and all
    datasets without a label source — stay `unlabeled_category`. When a fine
    column is not supplied for a source it falls back to the integration column.

    Returns (int_key, fine_key) — the two obs columns written on `adata`.
    """
    int_labels = pd.Series(unlabeled_category, index=adata.obs_names, dtype=object)
    fine_labels = pd.Series(unlabeled_category, index=adata.obs_names, dtype=object)
    bc_norm = adata.obs["barcodes"].map(_norm_barcode) if "barcodes" in adata.obs else \
        pd.Series(adata.obs_names, index=adata.obs_names).map(_norm_barcode)

    # Fall back to the integration column where no fine column is given.
    fine_cols = list(label_columns_fine or [])
    fine_cols += list(label_columns[len(fine_cols):])

    for dataset, ref_path, col_int, col_fine in zip(
            label_datasets, label_refs, label_columns, fine_cols):
        mask_ds = adata.obs[dataset_key].astype(str) == dataset
        n_ds = int(mask_ds.sum())
        if n_ds == 0:
            print(f"  [labels] dataset '{dataset}': no cells in merged object, skipping")
            continue
        if not os.path.exists(ref_path):
            print(f"  [labels] dataset '{dataset}': reference not found "
                  f"({ref_path}); leaving {n_ds} cells as '{unlabeled_category}'")
            continue

        ref = read_h5ad(ref_path)
        ref_bc = pd.Series(ref.obs_names, index=ref.obs_names).map(_norm_barcode)
        ds_bc = bc_norm[mask_ds]
        _apply_ref_column(ref, col_int, ref_bc, ds_bc, int_labels, dataset,
                          unlabeled_category)
        _apply_ref_column(ref, col_fine, ref_bc, ds_bc, fine_labels, dataset,
                          unlabeled_category)

    for name, labels, key in [("integration", int_labels, int_key),
                              ("fine", fine_labels, fine_key)]:
        n_labeled = int((labels != unlabeled_category).sum())
        print(f"  [labels] {name} ({key}): {n_labeled}/{adata.n_obs} labeled; "
              f"{labels[labels != unlabeled_category].nunique()} cell types")
        if n_labeled == 0:
            raise ValueError(
                f"scANVI selected but no cells could be labeled for the {name} "
                "column. Check BatchCorrection.scanvi.label_sources paths / "
                "column names (label_column / label_column_fine)."
            )
        adata.obs[key] = pd.Categorical(labels)

    return int_key, fine_key


def initialize_parser():
    parser = argparse.ArgumentParser(description='Batch correction (Harmony / scVI / scANVI)')
    parser.add_argument('--input', type=str, required=True, help="h5ad")
    parser.add_argument('--output', type=str, required=True, help="h5ad")
    parser.add_argument('--markers', type=str, required=True, help="json")
    parser.add_argument('--plotdir', type=str, required=True)
    parser.add_argument('--is_filtered', action='store_true')
    parser.add_argument('--sample_key', type=str, required=True)
    parser.add_argument('--donor_key', type=str, required=True)
    parser.add_argument('--dataset_key', type=str, default=None)
    parser.add_argument('--method', type=str, default='harmony',
                        choices=['harmony', 'scvi', 'scanvi'],
                        help="Batch correction method (scanvi always trains scVI first)")
    # scVI / scANVI hyperparameters
    parser.add_argument('--n_hvg', type=int, default=2000)
    parser.add_argument('--n_latent', type=int, default=30)
    parser.add_argument('--scvi_max_epochs', type=int, default=None,
                        help="None -> scvi-tools auto heuristic")
    parser.add_argument('--scanvi_max_epochs', type=int, default=20)
    parser.add_argument('--unlabeled_category', type=str, default='Unknown')
    # Per-dataset label sources for scANVI (parallel lists). `label_columns` is
    # the coarse integration label (scANVI #1); `label_columns_fine` is the fine
    # annotation label (scANVI #2) — both read from the harmonized reference.
    parser.add_argument('--label_datasets', type=str, nargs='*', default=[])
    parser.add_argument('--label_refs', type=str, nargs='*', default=[])
    parser.add_argument('--label_columns', type=str, nargs='*', default=[])
    parser.add_argument('--label_columns_fine', type=str, nargs='*', default=[])
    return parser


# -----------------------------
# Main
# -----------------------------
def main():
    parser = initialize_parser()
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    markers = args.markers
    plotdir = args.plotdir
    donor_key = args.donor_key
    sample_key = args.sample_key
    dataset_key = args.dataset_key
    is_filtered = args.is_filtered
    method = args.method

    sc.settings.figdir = plotdir

    print("1. Load Data")
    start = timeit.default_timer()
    adata = read_h5ad(input_file)
    stop = timeit.default_timer()
    print(f"Data loaded in {round(stop-start,2)}s")

    # Use dataset of origin (if present) + donor + sample as the batch covariates.
    batch_keys = [k for k in [dataset_key, donor_key, sample_key]
                  if k is not None and k in adata.obs.columns]

    print(f"2. Batch correction ({method})")
    start = timeit.default_timer()
    extra_colors = []
    stage_reps = []          # (obsm_key, tag) embeddings to plot per stage (scanvi)
    stage_label_col = None   # reference labels to colour those stage plots by
    if method == "harmony":
        use_rep = run_harmony(adata, batch_keys)
    elif method == "scvi":
        scvi_model = run_scvi(adata, batch_keys, args.n_hvg, args.n_latent,
                              args.scvi_max_epochs)
        use_rep = "X_scVI"
    elif method == "scanvi":
        # Method A (double scANVI): train scVI once, then fine-tune two heads on
        # it — #1 on the coarse integration labels (its embedding X_scANVI drives
        # clustering) and #2 on the fine labels (its predictions C_scANVI are the
        # annotation). The two label columns come from the harmonized references.
        int_key, fine_key = build_scanvi_labels(
            adata, dataset_key, args.label_datasets, args.label_refs,
            args.label_columns, args.label_columns_fine, args.unlabeled_category)
        scvi_model = run_scvi(adata, batch_keys, args.n_hvg, args.n_latent,
                              args.scvi_max_epochs)
        # #1 integration -> keep embedding (X_scANVI); coarse predictions kept for QC.
        run_scanvi(adata, scvi_model, int_key, args.unlabeled_category,
                   args.scanvi_max_epochs, rep_key="X_scANVI",
                   pred_key="C_scANVI_integration")
        # #2 annotation -> keep fine predictions (C_scANVI); its embedding
        # (X_scANVI_fine) is kept only for the per-stage diagnostic plot below.
        run_scanvi(adata, scvi_model, fine_key, args.unlabeled_category,
                   args.scanvi_max_epochs, rep_key="X_scANVI_fine",
                   pred_key="C_scANVI")
        use_rep = "X_scANVI"
        extra_colors = [int_key, "C_scANVI"]
        stage_reps = [("X_scVI", "scvi"),
                      ("X_scANVI", "scanvi_int"),
                      ("X_scANVI_fine", "scanvi_fine")]
        stage_label_col = fine_key
    else:
        raise ValueError(f"Unknown method '{method}'")
    stop = timeit.default_timer()
    print(f"Batches corrected in {round(stop-start,2)}s (rep: {use_rep})")

    if stage_reps:
        print("2b. Per-stage embedding plots (coloured by reference labels, no propagation)")
        start = timeit.default_timer()
        for rep, tag in stage_reps:
            plot_stage_embedding(adata, rep, stage_label_col, tag,
                                 args.unlabeled_category)
        stop = timeit.default_timer()
        print(f"Stage plots done in {round(stop-start,2)}s")

    print("3. Clustering")
    start = timeit.default_timer()
    clustering(adata, sample_key, donor_key, dataset_key, is_filtered, use_rep, extra_colors,
               unlabeled_category=args.unlabeled_category)
    stop = timeit.default_timer()
    print(f"Clustered in {round(stop-start,2)}s")

    print("4. Plot dotplots for Annotation")
    if markers != '':
        start = timeit.default_timer()
        adata = prepare_annotation(adata, markers)
        stop = timeit.default_timer()
        print(f"Dotplots plotted in in {round(stop-start,2)}s")
    else:
        print("No marker genes file provided")

    print("5. Write data")
    start = timeit.default_timer()
    # anndata refuses to write a frame whose index *name* also matches a column
    # (see MergeSamplesAnnData): the var index is named 'gene_symbols' while a
    # 'gene_symbols' column of un-suffixed symbols also exists. Drop the nominal
    # index name so the write succeeds.
    for frame in (adata.var, adata.obs):
        if frame.index.name in frame.columns:
            frame.index.name = None
    adata.write_h5ad(output_file)
    stop = timeit.default_timer()
    print(f"Data written in {round(stop-start,2)}s")


if __name__ == "__main__":
    start_total = timeit.default_timer()
    main()
    stop_total = timeit.default_timer()
    print(f"Total runtime: {round(stop_total-start_total,2)}s")
