#!/usr/bin/env python
"""
Genome-wide overlay of the SNP-array CNVs vs inferCNV, as three stacked bands:

    reference   SNP-array segments        (gain = red, loss = blue)
    query       inferCNV HMM regions      (gain = red, loss = blue)
    match       reference decomposed into  matched (green, same-direction query
                overlap) / contradicted (purple, opposite query) / missed (grey,
                no query call)

Chromosomes 1-22 are laid end-to-end along x (add X/Y with --include-xy). Both
inputs must be on hg38 (lift the array first with liftover_snp_array.py); inferCNV
is already hg38. Direction rules and parsing are shared with compare_cnv_overlap.py.

In analysis_mode "subclusters" the query has one cell_group_name per clone: pass
--group to restrict to one, or use plot_clones.py for the per-clone + summary set.

This module also exposes its drawing helpers (build_offsets, draw_dir_band,
draw_match_band, decompose_match, decorate_genome_axis, legend_handles,
render_overlay, colours) so plot_clones.py renders identical bands.

Usage:
    python plot_cnv_overlay.py \
        -r cellline.hg38.bed \
        -q results/BMO-SKNBE2c/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat \
        -o results/BMO-SKNBE2c/snp_vs_infercnv.overlay.png \
        --title BMO-SKNBE2c
"""

import os
import sys
import argparse

# reuse the exact parsing / direction logic from the overlap script (sibling file)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_cnv_overlap import parse_reference, parse_query, neutral_state_for  # noqa: E402


# GRCh38 primary-assembly chromosome sizes (bp).
HG38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}

C_GAIN = "#d6604d"      # red  - gain / amplification
C_LOSS = "#4393c3"      # blue - loss / deletion
C_MATCH = "#1a9850"     # green  - reference bp with same-direction inferCNV call
C_CONTRA = "#762a83"    # purple - reference bp with OPPOSITE inferCNV call
C_MISS = "#d9d9d9"      # grey   - reference event bp with no inferCNV call
C_SEP = "#bbbbbb"       # chromosome separators

ROW_Y = {"reference": 2.0, "query": 1.0, "match": 0.0}
ROW_H = 0.72


# --------------------------------------------------------------------------- #
# geometry / interval helpers
# --------------------------------------------------------------------------- #
def intersect(start, end, merged):
    """Clipped sub-intervals of [start,end) covered by sorted merged intervals."""
    out = []
    for s, e in merged:
        if e <= start:
            continue
        if s >= end:
            break
        out.append((max(start, s), min(end, e)))
    return out


def build_offsets(chroms):
    """Cumulative x-offset per chromosome, in the given order -> (offsets, total)."""
    off, cum = {}, 0
    for c in chroms:
        off[c] = cum
        cum += HG38[c]
    return off, cum


def bars(items, offsets):
    """(chrom,start,end) -> broken_barh xranges [(gx, width), ...], skipping
    chromosomes not on the plotted layout; returns (xranges, n_skipped)."""
    xr, skipped = [], 0
    for chrom, start, end in items:
        if chrom not in offsets:
            skipped += 1
            continue
        xr.append((offsets[chrom] + start, max(1, end - start)))
    return xr, skipped


def events_dir_items(events, direction):
    return [(e["chrom"], e["start"], e["end"]) for e in events if e["dir"] == direction]


def query_dir_items(query_by_dir, direction):
    return [(c, s, e) for c, ivs in query_by_dir[direction].items() for s, e in ivs]


def decompose_match(events, query_by_dir):
    """Split each reference event into missed / matched / contradicted bars and
    tally the bp. Returns a dict with bar lists, totals and pm/pc percentages."""
    opp = {"gain": "loss", "loss": "gain"}
    miss, match, contra = [], [], []
    tot_ref = tot_match = tot_contra = 0
    for e in events:
        chrom, s, t, d = e["chrom"], e["start"], e["end"], e["dir"]
        miss.append((chrom, s, t))                       # grey base = whole event
        same = intersect(s, t, query_by_dir[d].get(chrom, []))
        con = intersect(s, t, query_by_dir[opp[d]].get(chrom, []))
        match += [(chrom, a, b) for a, b in same]
        contra += [(chrom, a, b) for a, b in con]
        tot_ref += t - s
        tot_match += sum(b - a for a, b in same)
        tot_contra += sum(b - a for a, b in con)
    pm = 100.0 * tot_match / tot_ref if tot_ref else float("nan")
    pc = 100.0 * tot_contra / tot_ref if tot_ref else float("nan")
    return {"miss": miss, "match": match, "contra": contra,
            "tot_ref": tot_ref, "tot_match": tot_match, "tot_contra": tot_contra,
            "pm": pm, "pc": pc}


# --------------------------------------------------------------------------- #
# drawing helpers (shared with plot_clones.py)
# --------------------------------------------------------------------------- #
def _draw(ax, items, offsets, y, h, color):
    xr, sk = bars(items, offsets)
    if xr:
        ax.broken_barh(xr, (y, h), facecolors=color, edgecolors="none")
    return sk


def draw_dir_band(ax, gain_items, loss_items, offsets, y, h=ROW_H):
    """A gain/loss band (reference segments or query intervals)."""
    return (_draw(ax, gain_items, offsets, y, h, C_GAIN)
            + _draw(ax, loss_items, offsets, y, h, C_LOSS))


def draw_match_band(ax, events, query_by_dir, offsets, y, h=ROW_H):
    """A match band: grey base, then matched (green) and contradicted (purple)
    overlaid. Returns decompose_match() dict with an added 'skipped' count."""
    dec = decompose_match(events, query_by_dir)
    sk = _draw(ax, dec["miss"], offsets, y, h, C_MISS)
    sk += _draw(ax, dec["match"], offsets, y, h, C_MATCH)
    sk += _draw(ax, dec["contra"], offsets, y, h, C_CONTRA)
    dec["skipped"] = sk
    return dec


def decorate_genome_axis(ax, chroms, offsets, genome_len, fontsize=8):
    """Chromosome separators, centred x tick labels and x limits."""
    for c in chroms:
        ax.axvline(offsets[c], color=C_SEP, lw=0.6, zorder=0)
    ax.axvline(genome_len, color=C_SEP, lw=0.6, zorder=0)
    ax.set_xticks([offsets[c] + HG38[c] / 2 for c in chroms])
    ax.set_xticklabels([c.replace("chr", "") for c in chroms], fontsize=fontsize)
    ax.set_xlim(0, genome_len)


def legend_handles(match_only=False):
    from matplotlib.patches import Patch
    match = [Patch(facecolor=C_MATCH, label="matched (same dir)"),
             Patch(facecolor=C_CONTRA, label="contradicted (opp dir)"),
             Patch(facecolor=C_MISS, label="missed (no call)")]
    if match_only:
        return match
    return [Patch(facecolor=C_GAIN, label="gain"),
            Patch(facecolor=C_LOSS, label="loss")] + match


def genome_chroms(include_xy=False):
    return [f"chr{i}" for i in range(1, 23)] + (["chrX", "chrY"] if include_xy else [])


def import_pyplot():
    """Import matplotlib (Agg backend), with an actionable message if it's absent."""
    try:
        import matplotlib
    except ImportError:
        sys.exit(
            "[plot] matplotlib is not installed in this Python. Fix with one of:\n"
            "         pip install matplotlib\n"
            "         conda install -c conda-forge matplotlib-base\n"
            "       or run in an env that already has it (the pipeline env\n"
            "       inferCNV/snp_array/envs/snp_array.yaml, or your scverse env).")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def render_overlay(out, events, query_by_dir, title="", chroms=None,
                   width=18.0, height=4.2, dpi=150, ref_label="SNP array"):
    """Draw + save the 3-band overlay for one (reference, query) pair. `query_by_dir`
    is a {dir:{chrom:[merged]}} (e.g. a single clone). `ref_label` names the reference
    (e.g. 'SNP array' or another sample's inferCNV). Returns the decompose dict (with
    'skipped'). Shared by the CLI and plot_clones.py."""
    plt = import_pyplot()

    if chroms is None:
        chroms = genome_chroms(False)
    offsets, genome_len = build_offsets(chroms)

    fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)
    sk = draw_dir_band(ax, events_dir_items(events, "gain"),
                       events_dir_items(events, "loss"), offsets, ROW_Y["reference"])
    sk += draw_dir_band(ax, query_dir_items(query_by_dir, "gain"),
                        query_dir_items(query_by_dir, "loss"), offsets, ROW_Y["query"])
    dec = draw_match_band(ax, events, query_by_dir, offsets, ROW_Y["match"])
    sk += dec["skipped"]

    decorate_genome_axis(ax, chroms, offsets, genome_len)
    ax.set_yticks([ROW_Y[r] + ROW_H / 2 for r in ("match", "query", "reference")])
    ax.set_yticklabels(["match", "query\n(inferCNV)", f"reference\n({ref_label})"], fontsize=9)
    ax.set_ylim(-0.25, ROW_Y["reference"] + ROW_H + 0.25)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("chromosome", fontsize=9)

    t = title + "  " if title else ""
    ax.set_title(f"{t}{ref_label} vs inferCNV  —  "
                 f"{dec['pm']:.0f}% of reference bp matched, {dec['pc']:.0f}% contradicted",
                 fontsize=11)
    ax.legend(handles=legend_handles(), ncol=5, fontsize=8, loc="lower center",
              bbox_to_anchor=(0.5, -0.32), frameon=False)

    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    dec["skipped_total"] = sk
    return dec


# --------------------------------------------------------------------------- #
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-r", "--reference", required=True, help="SNP-array CNV BED (hg38)")
    p.add_argument("-q", "--query", required=True, help="inferCNV pred_cnv_regions.dat")
    p.add_argument("-o", "--out", default="cnv_overlay.png", help="output image (.png/.pdf)")
    p.add_argument("--type-col", type=int, default=4,
                   help="1-based category column for event direction (default 4)")
    p.add_argument("--cn-col", type=int, default=13,
                   help="1-based copy-number column, recorded only (default 13)")
    p.add_argument("--group", default=None, help="restrict inferCNV to one cell_group_name")
    p.add_argument("--hmm-i", type=int, choices=(3, 6), default=None,
                   help="inferCNV HMM model (i3/i6); default: autodetect from filename")
    p.add_argument("--include-xy", action="store_true", help="also plot chrX/chrY")
    p.add_argument("--title", default="", help="plot title (e.g. the sample name)")
    p.add_argument("--width", type=float, default=18.0)
    p.add_argument("--height", type=float, default=4.2)
    p.add_argument("--dpi", type=int, default=150)
    args = p.parse_args(argv)

    events = parse_reference(args.reference, args.type_col, args.cn_col)
    neutral_state, hmm_i = neutral_state_for(args.query, args.hmm_i)
    print(f"[plot] inferCNV model: HMMi{hmm_i} (neutral = state {neutral_state})")
    query = parse_query(args.query, args.group, neutral_state)
    dec = render_overlay(args.out, events, query, title=args.title,
                         chroms=genome_chroms(args.include_xy),
                         width=args.width, height=args.height, dpi=args.dpi)

    if dec["skipped_total"]:
        print(f"[plot] note: skipped {dec['skipped_total']} segment(s) on chromosomes "
              f"not plotted ({'chr1-22' if not args.include_xy else 'chr1-22,X,Y'}).")
    print(f"[plot] reference bp matched: {dec['pm']:.1f}%  contradicted: {dec['pc']:.1f}%")
    print(f"[plot] wrote {args.out}")


if __name__ == "__main__":
    main()
