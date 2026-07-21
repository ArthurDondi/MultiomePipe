#!/usr/bin/env python
"""
Per-clone SNP-array vs inferCNV comparison for one sample run in analysis_mode
"subclusters", where each cell_group_name in pred_cnv_regions.dat is a clone.

Produces, taking the SNP array as reference:
  1. one 3-band overlay per clone  (reference / that clone / match)   [--no-overlays to skip]
  2. a SUMMARY plot stacking every clone's match band on shared genome axes, with a
     clone-size sidebar (cells per clone, from observation_groupings.txt)
  3. a per-clone table  <prefix>.clones.tsv

Clone sizes come from inferCNV's observation_groupings.txt: the "Dendrogram Group"
column is the subcluster id (e.g. neuron_s8) and matches the suffix of a
pred_cnv_regions cell_group_name (e.g. neuron.neuron_s8). Rows in the summary are
sorted largest-clone first, so a big well-matching clone stands out from tiny noisy
ones. Everything must be hg38 (lift the array first with liftover_snp_array.py).

Usage (one sample):
    python plot_clones.py \
        -r cellline_SKNBE2c.hg38.bed \
        -q results/BMO-SKNBE2c/17_HMM_predHMMi3.leiden.hmm_mode-subclusters.pred_cnv_regions.dat \
        -g results/BMO-SKNBE2c/infercnv.17_HMM_predHMMi3.leiden.hmm_mode-subclusters.observation_groupings.txt \
        -o results/BMO-SKNBE2c/snp_vs_infercnv \
        --title BMO-SKNBE2c

Loop over samples with a shell for-loop (see README).
"""

import os
import re
import sys
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_cnv_overlap import (  # noqa: E402
    parse_reference, parse_query_by_group, parse_groupings, clone_size,
    compare, summarize, neutral_state_for)
from plot_cnv_overlay import (  # noqa: E402
    ROW_H, C_GAIN, C_LOSS, build_offsets, genome_chroms, events_dir_items,
    draw_dir_band, draw_match_band, decorate_genome_axis, legend_handles,
    render_overlay)

C_SIZE = "#9e9ac8"     # clone-size sidebar bars


def safe(name):
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def short(name):
    """Drop the annotation prefix: 'neuron.neuron_s8' -> 'neuron_s8'."""
    return name.split(".", 1)[1] if "." in name else name


def per_clone_metrics(events, by_group, sizes, min_overlap):
    """One record per clone: size + event/bp match metrics (SNP array as ref)."""
    recs = []
    for clone, qbd in by_group.items():
        summ = summarize(compare(events, qbd, min_overlap), min_overlap)["all"]
        recs.append({
            "clone": clone,
            "n_cells": clone_size(clone, sizes),
            "n_events": summ["n_events"],
            "n_events_matched": summ["n_matched"],
            "pct_events_matched": summ["pct_events_matched"],
            "pct_bp_matched": summ["pct_bp_matched"],
            "pct_bp_contradicted": summ["pct_bp_contradicted"],
        })
    return recs


def sort_key(rec, has_sizes):
    # largest clone first when sizes are known; else best bp-match first
    if has_sizes:
        return (-(rec["n_cells"] if rec["n_cells"] is not None else -1),
                -_num(rec["pct_bp_matched"]))
    return (-_num(rec["pct_bp_matched"]),)


def _num(x):
    return -1.0 if x != x else x        # NaN -> -1 for sorting


def write_table(recs, path):
    with open(path, "w") as fh:
        fh.write("clone\tn_cells\tn_events\tn_events_matched\tpct_events_matched\t"
                 "pct_bp_matched\tpct_bp_contradicted\n")
        for r in recs:
            fh.write(f"{r['clone']}\t{'' if r['n_cells'] is None else r['n_cells']}\t"
                     f"{r['n_events']}\t{r['n_events_matched']}\t"
                     f"{r['pct_events_matched']:.2f}\t{r['pct_bp_matched']:.2f}\t"
                     f"{r['pct_bp_contradicted']:.2f}\n")
    print(f"[clones] wrote {path}")


def summary_plot(events, by_group, recs, out, chroms, title, width, dpi):
    """Stacked match band per clone (sorted) + a clone-size sidebar."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    offsets, genome_len = build_offsets(chroms)
    has_sizes = any(r["n_cells"] is not None for r in recs)
    total_cells = sum(r["n_cells"] for r in recs if r["n_cells"] is not None)

    # rows top -> bottom: reference context row, then clones (already sorted)
    order = ["__ref__"] + [r["clone"] for r in recs]
    n = len(order)
    y_of = {name: (n - 1 - i) for i, name in enumerate(order)}   # top row highest y
    rec_by_clone = {r["clone"]: r for r in recs}

    fig_h = max(3.0, 0.34 * n + 1.4)
    if has_sizes:
        fig, (ax, ax2) = plt.subplots(
            1, 2, figsize=(width + 2.2, fig_h), sharey=True,
            gridspec_kw={"width_ratios": [width, 2.0], "wspace": 0.03},
            constrained_layout=True)
    else:
        fig, ax = plt.subplots(figsize=(width, fig_h), constrained_layout=True)
        ax2 = None

    # reference context row
    draw_dir_band(ax, events_dir_items(events, "gain"),
                  events_dir_items(events, "loss"), offsets, y_of["__ref__"])
    # one match band per clone
    for clone in order[1:]:
        draw_match_band(ax, events, by_group[clone], offsets, y_of[clone])

    decorate_genome_axis(ax, chroms, offsets, genome_len)
    yt, yl = [], []
    for name in order:
        yt.append(y_of[name] + ROW_H / 2)
        if name == "__ref__":
            yl.append("reference (SNP array)")
        else:
            r = rec_by_clone[name]
            pm = r["pct_bp_matched"]
            yl.append(f"{short(name)}  {0 if pm != pm else round(pm)}%")
    ax.set_yticks(yt)
    ax.set_yticklabels(yl, fontsize=7)
    ax.set_ylim(-0.4, n - 1 + ROW_H + 0.4)
    ax.tick_params(axis="y", length=0)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.set_xlabel("chromosome", fontsize=9)

    t = (title + "  ") if title else ""
    ax.set_title(f"{t}per-clone match to SNP array  "
                 f"(row label = % of reference bp matched)", fontsize=11)
    ax.legend(handles=legend_handles(), ncol=5, fontsize=8, loc="lower center",
              bbox_to_anchor=(0.5, -0.9 / fig_h - 0.02), frameon=False)

    if ax2 is not None:
        centers = [y_of[r["clone"]] + ROW_H / 2 for r in recs if r["n_cells"] is not None]
        vals = [r["n_cells"] for r in recs if r["n_cells"] is not None]
        ax2.barh(centers, vals, height=ROW_H, color=C_SIZE)
        vmax = max(vals) if vals else 1
        for c, v in zip(centers, vals):
            ax2.text(v + vmax * 0.02, c, str(v), va="center", fontsize=6)
        ax2.set_xlim(0, vmax * 1.18)
        ax2.set_xlabel("cells", fontsize=8)
        ax2.tick_params(axis="x", labelsize=7)
        for spine in ("top", "right", "left"):
            ax2.spines[spine].set_visible(False)
        ax2.set_title(f"clone size\n(total {total_cells})", fontsize=8)

    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    print(f"[clones] wrote {out}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-r", "--reference", required=True, help="SNP-array CNV BED (hg38)")
    p.add_argument("-q", "--query", required=True,
                   help="inferCNV subclusters pred_cnv_regions.dat")
    p.add_argument("-g", "--groupings", default=None,
                   help="inferCNV observation_groupings.txt (for clone sizes)")
    p.add_argument("-o", "--out-prefix", required=True,
                   help="output path prefix (writes <prefix>.clones_summary.png, "
                        ".clones.tsv, and <prefix>.<clone>.overlay.png)")
    p.add_argument("--cn-col", type=int, default=13)
    p.add_argument("--neutral-cn", type=float, default=2.0)
    p.add_argument("--hmm-i", type=int, choices=(3, 6), default=None,
                   help="inferCNV HMM model (i3/i6); default: autodetect from filename")
    p.add_argument("--min-overlap", type=float, default=0.5)
    p.add_argument("--min-cells", type=int, default=0,
                   help="drop clones with fewer than this many cells (default 0 = keep all)")
    p.add_argument("--include-xy", action="store_true")
    p.add_argument("--title", default="", help="sample name for titles")
    p.add_argument("--no-overlays", action="store_true",
                   help="skip the per-clone overlay PNGs (summary + table only)")
    p.add_argument("--width", type=float, default=18.0)
    p.add_argument("--dpi", type=int, default=150)
    args = p.parse_args(argv)

    chroms = genome_chroms(args.include_xy)
    events = parse_reference(args.reference, args.cn_col, args.neutral_cn)
    if not events:
        sys.exit("[clones] no gain/loss events in the reference — nothing to compare.")
    neutral_state, hmm_i = neutral_state_for(args.query, args.hmm_i)
    print(f"[clones] inferCNV model: HMMi{hmm_i} (neutral = state {neutral_state})")
    by_group = parse_query_by_group(args.query, neutral_state=neutral_state)
    sizes = parse_groupings(args.groupings) if args.groupings else {}

    recs = per_clone_metrics(events, by_group, sizes, args.min_overlap)
    has_sizes = any(r["n_cells"] is not None for r in recs)

    if args.min_cells > 0:
        keep = [r for r in recs if (r["n_cells"] or 0) >= args.min_cells]
        dropped = len(recs) - len(keep)
        if dropped:
            print(f"[clones] dropping {dropped} clone(s) with < {args.min_cells} cells")
        recs = keep
    if not recs:
        sys.exit("[clones] no clones left to plot (check --min-cells / inputs).")

    recs.sort(key=lambda r: sort_key(r, has_sizes))
    by_group = {r["clone"]: by_group[r["clone"]] for r in recs}   # same order

    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)
    write_table(recs, args.out_prefix + ".clones.tsv")

    if not args.no_overlays:
        for i, r in enumerate(recs, 1):
            clone = r["clone"]
            size = "" if r["n_cells"] is None else f", n={r['n_cells']}"
            ttl = f"{args.title + ' · ' if args.title else ''}{clone}{size}"
            out = f"{args.out_prefix}.{safe(clone)}.overlay.png"
            render_overlay(out, events, by_group[clone], title=ttl, chroms=chroms,
                           width=args.width, dpi=args.dpi)
            print(f"[clones]   overlay {i}/{len(recs)}: {out}")

    summary_plot(events, by_group, recs, args.out_prefix + ".clones_summary.png",
                 chroms, args.title, args.width, args.dpi)
    print(f"[clones] done: {len(recs)} clone(s).")


if __name__ == "__main__":
    main()
