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
from compare_cnv_overlap import parse_reference, parse_query, norm_chrom  # noqa: E402


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
    """Cumulative x-offset per chromosome, in the given order."""
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


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-r", "--reference", required=True, help="SNP-array CNV BED (hg38)")
    p.add_argument("-q", "--query", required=True, help="inferCNV pred_cnv_regions.dat")
    p.add_argument("-o", "--out", default="cnv_overlay.png", help="output image (.png/.pdf)")
    p.add_argument("--cn-col", type=int, default=13, help="1-based CN column (default 13)")
    p.add_argument("--neutral-cn", type=float, default=2.0)
    p.add_argument("--group", default=None, help="restrict inferCNV to one cell_group_name")
    p.add_argument("--include-xy", action="store_true", help="also plot chrX/chrY")
    p.add_argument("--title", default="", help="plot title (e.g. the sample name)")
    p.add_argument("--width", type=float, default=18.0)
    p.add_argument("--height", type=float, default=4.2)
    p.add_argument("--dpi", type=int, default=150)
    args = p.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    chroms = [f"chr{i}" for i in range(1, 23)]
    if args.include_xy:
        chroms += ["chrX", "chrY"]
    offsets, genome_len = build_offsets(chroms)

    events = parse_reference(args.reference, args.cn_col, args.neutral_cn)
    query = parse_query(args.query, args.group)

    # ---- reference band ----
    ref_gain = [(e["chrom"], e["start"], e["end"]) for e in events if e["dir"] == "gain"]
    ref_loss = [(e["chrom"], e["start"], e["end"]) for e in events if e["dir"] == "loss"]

    # ---- query band (merged intervals per direction) ----
    q_gain = [(c, s, e) for c, ivs in query["gain"].items() for s, e in ivs]
    q_loss = [(c, s, e) for c, ivs in query["loss"].items() for s, e in ivs]

    # ---- match band: decompose each reference event ----
    opp = {"gain": "loss", "loss": "gain"}
    miss_bars, match_bars, contra_bars = [], [], []
    tot_ref = tot_match = tot_contra = 0
    for e in events:
        chrom, s, t, d = e["chrom"], e["start"], e["end"], e["dir"]
        miss_bars.append((chrom, s, t))                       # grey base = whole event
        same = intersect(s, t, query[d].get(chrom, []))
        contra = intersect(s, t, query[opp[d]].get(chrom, []))
        match_bars += [(chrom, a, b) for a, b in same]
        contra_bars += [(chrom, a, b) for a, b in contra]
        tot_ref += t - s
        tot_match += sum(b - a for a, b in same)
        tot_contra += sum(b - a for a, b in contra)

    # ---- draw ----
    fig, ax = plt.subplots(figsize=(args.width, args.height), constrained_layout=True)
    skipped_total = 0

    def draw(items, row, color):
        nonlocal skipped_total
        xr, sk = bars(items, offsets)
        skipped_total += sk
        if xr:
            ax.broken_barh(xr, (ROW_Y[row], ROW_H), facecolors=color,
                           edgecolors="none")

    draw(ref_gain, "reference", C_GAIN)
    draw(ref_loss, "reference", C_LOSS)
    draw(q_gain, "query", C_GAIN)
    draw(q_loss, "query", C_LOSS)
    draw(miss_bars, "match", C_MISS)       # grey base first ...
    draw(match_bars, "match", C_MATCH)     # ... then overlay matched (green) ...
    draw(contra_bars, "match", C_CONTRA)   # ... and contradicted (purple) on top

    # chromosome separators + centered labels
    for c in chroms:
        x0 = offsets[c]
        ax.axvline(x0, color=C_SEP, lw=0.6, zorder=0)
    ax.axvline(genome_len, color=C_SEP, lw=0.6, zorder=0)
    ax.set_xticks([offsets[c] + HG38[c] / 2 for c in chroms])
    ax.set_xticklabels([c.replace("chr", "") for c in chroms], fontsize=8)
    ax.set_xlim(0, genome_len)

    ax.set_yticks([ROW_Y[r] + ROW_H / 2 for r in ("match", "query", "reference")])
    ax.set_yticklabels(["match", "query\n(inferCNV)", "reference\n(SNP array)"], fontsize=9)
    ax.set_ylim(-0.25, ROW_Y["reference"] + ROW_H + 0.25)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("chromosome", fontsize=9)

    pm = 100.0 * tot_match / tot_ref if tot_ref else float("nan")
    pc = 100.0 * tot_contra / tot_ref if tot_ref else float("nan")
    title = args.title + "  " if args.title else ""
    ax.set_title(f"{title}SNP array vs inferCNV  —  "
                 f"{pm:.0f}% of reference bp matched, {pc:.0f}% contradicted",
                 fontsize=11)

    legend = [Patch(facecolor=C_GAIN, label="gain"),
              Patch(facecolor=C_LOSS, label="loss"),
              Patch(facecolor=C_MATCH, label="matched (same dir)"),
              Patch(facecolor=C_CONTRA, label="contradicted (opp dir)"),
              Patch(facecolor=C_MISS, label="missed (no call)")]
    ax.legend(handles=legend, ncol=5, fontsize=8, loc="lower center",
              bbox_to_anchor=(0.5, -0.32), frameon=False)

    if skipped_total:
        print(f"[plot] note: skipped {skipped_total} segment(s) on chromosomes not "
              f"plotted ({'chr1-22' if not args.include_xy else 'chr1-22,X,Y'}).")

    fig.savefig(args.out, dpi=args.dpi)
    print(f"[plot] reference bp matched: {pm:.1f}%  contradicted: {pc:.1f}%")
    print(f"[plot] wrote {args.out}")


if __name__ == "__main__":
    main()
