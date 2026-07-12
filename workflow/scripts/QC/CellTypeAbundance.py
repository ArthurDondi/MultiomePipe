import os
import json
import argparse
import timeit
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Cell-type list ("frenz_ctypes")
# -----------------------------------------------------------------------------
def load_celltypes(celltypes_json=None, celltypes_inline=None,
                   reference=None, reference_celltype_key="author_cell_type"):
    """Return the premade list of cell types to restrict the comparison to.

    Resolved in priority order so the same script works interactively and on
    the CLI:
      1. ``celltypes_inline`` — a Python list passed directly (interactive use).
      2. ``celltypes_json``   — path to a JSON file holding either a plain
         array ``["A", "B", ...]`` or an object ``{"frenz_ctypes": [...]}``.
      3. ``reference``        — a Frenz reference h5ad; the list is derived from
         the categories of ``reference_celltype_key`` (literally "the cell types
         present in the Frenz reference").
    """
    if celltypes_inline:
        return list(celltypes_inline)

    if celltypes_json:
        with open(celltypes_json) as f:
            data = json.load(f)
        if isinstance(data, dict):
            # Accept {"frenz_ctypes": [...]} or the first list-valued entry.
            data = data.get(
                "frenz_ctypes",
                next((v for v in data.values() if isinstance(v, list)), None),
            )
        if not isinstance(data, list):
            raise ValueError(
                f"{celltypes_json} must contain a JSON array of cell-type names "
                "(or an object with a 'frenz_ctypes' array)."
            )
        return list(data)

    if reference:
        ref = ad.io.read_h5ad(reference)
        if reference_celltype_key not in ref.obs.columns:
            raise ValueError(
                f"'{reference_celltype_key}' not in reference obs. "
                f"Available: {list(ref.obs.columns)}"
            )
        col = ref.obs[reference_celltype_key].astype("category")
        return col.cat.categories.astype(str).tolist()

    raise ValueError(
        "No cell-type list provided. Pass --celltypes (JSON), --celltype "
        "(inline, repeatable), or --reference to derive it from a Frenz h5ad."
    )


# -----------------------------------------------------------------------------
# Grouping
# -----------------------------------------------------------------------------
def assign_groups(adata, sample_key, group_specs):
    """Label each cell with a group based on a substring of its sample name.

    ``group_specs`` is an ordered list of ``(name, pattern)`` tuples; a cell is
    assigned to the first group whose ``pattern`` is a substring of its
    ``sample_key`` value. Cells matching no group are left as NA (dropped).
    Substring matching is exactly what "all samples containing 'Frenz'" asks
    for; 'BMOalone' only ever matches the BMOalone sample.
    """
    if sample_key not in adata.obs.columns:
        raise ValueError(
            f"'{sample_key}' not in obs. Available: {list(adata.obs.columns)}"
        )
    samples = adata.obs[sample_key].astype(str)
    group = pd.Series(pd.NA, index=adata.obs_names, dtype=object)
    for name, pattern in group_specs:
        mask = samples.str.contains(pattern, regex=False) & group.isna()
        group[mask] = name
    return group


# -----------------------------------------------------------------------------
# Composition table
# -----------------------------------------------------------------------------
def composition_table(adata, group_key, celltype_key, frenz_ctypes, group_order):
    """cell-type x group count matrix, rows/cols ordered, restricted to list."""
    if celltype_key not in adata.obs.columns:
        raise ValueError(
            f"'{celltype_key}' not in obs. Available: {list(adata.obs.columns)}"
        )
    df = adata.obs[[group_key, celltype_key]].copy()
    df = df[df[group_key].notna()]
    df[celltype_key] = df[celltype_key].astype(str)

    # Restrict to the premade cell-type list.
    present = [c for c in frenz_ctypes if c in set(df[celltype_key])]
    missing = [c for c in frenz_ctypes if c not in present]
    dropped_ct = sorted(set(df[celltype_key]) - set(frenz_ctypes))
    df = df[df[celltype_key].isin(frenz_ctypes)]

    counts = (
        df.groupby([group_key, celltype_key]).size().unstack(fill_value=0)
        .reindex(index=[g for g in group_order if g in set(df[group_key])],
                 columns=[c for c in frenz_ctypes if c in present],
                 fill_value=0)
    )
    return counts, present, missing, dropped_ct


# -----------------------------------------------------------------------------
# Colours
# -----------------------------------------------------------------------------
def build_color_map(adata, celltype_key, celltypes):
    """Reuse the pipeline's per-cell-type colours when available, else palette."""
    base = {}
    uns_key = f"{celltype_key}_colors"
    col = adata.obs.get(celltype_key)
    if (col is not None and hasattr(col, "cat")
            and uns_key in adata.uns
            and len(adata.uns[uns_key]) == len(col.cat.categories)):
        base = dict(zip(col.cat.categories.astype(str), list(adata.uns[uns_key])))

    palette = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    for i, c in enumerate([c for c in celltypes if c not in base]):
        base[c] = matplotlib.colors.to_hex(palette[i % len(palette)])
    return {c: base[c] for c in celltypes}


# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
def _stacked_bar(mat, color_map, outpath, ylabel, title, ylim=None):
    groups = list(mat.index)
    x = np.arange(len(groups))
    fig, ax = plt.subplots(figsize=(max(3.5, 1.7 * len(groups) + 2), 6))
    bottom = np.zeros(len(groups))
    for ct in mat.columns:
        vals = mat[ct].to_numpy(dtype=float)
        ax.bar(x, vals, bottom=bottom, width=0.7, label=ct,
               color=color_map[ct], edgecolor="white", linewidth=0.3)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=0, ha="center")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.legend(title="Cell type", bbox_to_anchor=(1.02, 1), loc="upper left",
              frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_celltype_abundance(adata, frenz_ctypes, plotdir=".",
                            celltype_key="author_cell_type", sample_key="sample",
                            group_specs=(("Frenz", "Frenz"), ("BMOalone", "BMOalone")),
                            prefix="celltype_abundance_frenz_vs_bmoalone",
                            output_csv=None):
    """Stacked barplots of cell-type composition for the requested groups.

    Produces a proportion (100%-stacked) plot and a raw-count plot, restricted
    to ``frenz_ctypes``. Returns the count table. Importable for interactive use:

        plot_celltype_abundance(adata, frenz_ctypes)
    """
    os.makedirs(plotdir, exist_ok=True)
    group_specs = [tuple(g) for g in group_specs]
    group_order = [name for name, _ in group_specs]

    adata.obs["group"] = assign_groups(adata, sample_key, group_specs)
    counts, present, missing, dropped_ct = composition_table(
        adata, "group", celltype_key, frenz_ctypes, group_order
    )

    if counts.empty or counts.to_numpy().sum() == 0:
        raise ValueError(
            "No cells matched the requested groups AND cell-type list. "
            f"Groups tried: {group_specs}; cell-type key: '{celltype_key}'."
        )

    # ---- report -------------------------------------------------------------
    n_per_group = counts.sum(axis=1)
    print("Cells per group (within frenz_ctypes):")
    for g in counts.index:
        print(f"  {g}: {int(n_per_group[g])}")
    if missing:
        print(f"frenz_ctypes not observed in the selected cells ({len(missing)}): {missing}")
    if dropped_ct:
        print(f"Observed cell types NOT in frenz_ctypes (excluded): {dropped_ct}")

    color_map = build_color_map(adata, celltype_key, list(counts.columns))

    # ---- proportions (the honest cross-group comparison) --------------------
    props = counts.div(counts.sum(axis=1), axis=0)
    labelled = props.copy()
    labelled.index = [f"{g}\n(n={int(n_per_group[g])})" for g in props.index]
    _stacked_bar(
        labelled, color_map,
        os.path.join(plotdir, f"{prefix}_proportion.png"),
        ylabel="Fraction of cells", ylim=(0, 1),
        title="Cell-type composition (Frenz cell types only)",
    )

    # ---- raw counts ---------------------------------------------------------
    labelled_counts = counts.copy()
    labelled_counts.index = [f"{g}\n(n={int(n_per_group[g])})" for g in counts.index]
    _stacked_bar(
        labelled_counts, color_map,
        os.path.join(plotdir, f"{prefix}_counts.png"),
        ylabel="Number of cells", title="Cell-type counts (Frenz cell types only)",
    )

    # ---- table --------------------------------------------------------------
    out = counts.T.copy()
    out.columns = [f"{g}_count" for g in out.columns]
    for g in counts.index:
        out[f"{g}_fraction"] = (counts.loc[g] / n_per_group[g]).values
    out.index.name = celltype_key
    if output_csv:
        out.to_csv(output_csv)
        print(f"Wrote composition table: {output_csv}")

    print(f"Wrote plots to {plotdir}/{prefix}_proportion.png and {prefix}_counts.png")
    return counts


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def _parse_group_specs(items):
    specs = []
    for it in items:
        if "=" not in it:
            raise ValueError(f"--groups entry '{it}' must be name=pattern")
        name, pattern = it.split("=", 1)
        specs.append((name, pattern))
    return specs


def initialize_parser():
    parser = argparse.ArgumentParser(
        description="Stacked barplots of cell-type abundance between sample "
                    "groups, restricted to a premade cell-type list (frenz_ctypes)."
    )
    parser.add_argument("--input", required=True,
                        help="Query h5ad (e.g. merged.label_transferred.h5ad)")
    parser.add_argument("--plotdir", required=True, help="Output directory for plots")
    parser.add_argument("--output_csv", default=None, help="Optional composition CSV")

    # frenz_ctypes (choose one source)
    parser.add_argument("--celltypes", default=None,
                        help="JSON file with the frenz_ctypes list")
    parser.add_argument("--celltype", action="append", default=None,
                        help="Inline cell-type name; repeat for several")
    parser.add_argument("--reference", default=None,
                        help="Frenz reference h5ad to derive the list from")
    parser.add_argument("--reference_celltype_key", default="author_cell_type",
                        help="Reference obs column for --reference (default: author_cell_type)")

    parser.add_argument("--celltype_key", default="author_cell_type",
                        help="obs column with cell types to plot (default: author_cell_type)")
    parser.add_argument("--sample_key", default="sample",
                        help="obs column with sample names (default: sample)")
    parser.add_argument("--groups", nargs="+", default=["Frenz=Frenz", "BMOalone=BMOalone"],
                        help="Ordered name=substring specs (default: Frenz=Frenz BMOalone=BMOalone)")
    return parser


def main():
    args = initialize_parser().parse_args()

    frenz_ctypes = load_celltypes(
        celltypes_json=args.celltypes,
        celltypes_inline=args.celltype,
        reference=args.reference,
        reference_celltype_key=args.reference_celltype_key,
    )
    print(f"frenz_ctypes ({len(frenz_ctypes)}): {frenz_ctypes}")

    print("Load query data")
    adata = ad.io.read_h5ad(args.input)

    plot_celltype_abundance(
        adata,
        frenz_ctypes,
        plotdir=args.plotdir,
        celltype_key=args.celltype_key,
        sample_key=args.sample_key,
        group_specs=_parse_group_specs(args.groups),
        output_csv=args.output_csv,
    )


if __name__ == "__main__":
    start = timeit.default_timer()
    main()
    print(f"Total runtime: {round(timeit.default_timer() - start, 2)}s")
