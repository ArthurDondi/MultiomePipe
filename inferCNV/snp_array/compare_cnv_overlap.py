#!/usr/bin/env python
"""
Overlap a cell-line SNP-array CNV BED (hg38) against an inferCNV HMM prediction
file, taking the SNP array as the reference / ground truth.

    reference = SNP array  (pred_cnv_regions-style truth for the cell line)
    query     = inferCNV   (17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat)

Both files must be on the SAME genome build. The array is hg19, so lift it first
with liftover_snp_array.py; inferCNV is already hg38.

Direction of an event (only the direction is compared, never absolute copies):
  reference (SNP array): copy-number column, default col 13 (1-based), values ~1-40
        cn > 2  -> gain      cn < 2  -> loss      cn == 2 -> neutral (ignored)
  query     (inferCNV):  HMM `state` column, 1..6 with 3 = diploid baseline
        state > 3 -> gain    state < 3 -> loss    state == 3 -> neutral (never emitted)

A reference event is a MATCH when same-direction query calls cover at least
--min-overlap of its length (default 0.50). Only same-chromosome, same-direction
base pairs count. Overlapping query intervals are merged per direction first, so
base pairs are never double-counted (this also unions several inferCNV cell
groups / subclusters cleanly).

Two headline numbers, both with the SNP array as reference:
  1. % events matching = matched reference events / all reference events
  2. % bp matching     = same-direction overlapping bp / all reference-event bp
     -> the denominator includes reference events that do NOT clear the per-event
        match threshold, i.e. it is a pure base-pair directional concordance.
Both are reported overall and split into gains / losses. Opposite-direction
("contradicted") bp are also reported as QC.

A per-event table (--out-prefix -> *.per_event.tsv) lists every reference event
with its same/opposite overlap so misses and disagreements can be inspected.

Usage:
    python compare_cnv_overlap.py \
        -r cellline.hg38.bed \
        -q results/BMO-SKNBE2c/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat \
        --out-prefix results/BMO-SKNBE2c/snp_vs_infercnv
"""

import sys
import argparse
from collections import defaultdict


# --------------------------------------------------------------------------- #
# parsing
# --------------------------------------------------------------------------- #
def norm_chrom(c):
    return c if c.startswith("chr") else "chr" + c


def cn_to_dir(cn, neutral_cn):
    if cn > neutral_cn:
        return "gain"
    if cn < neutral_cn:
        return "loss"
    return "neutral"


def state_to_dir(state):
    # inferCNV i6 HMM: 1,2 = loss (0x, 0.5x); 3 = diploid; 4,5,6 = gain (1.5x,2x,3x)
    if state > 3:
        return "gain"
    if state < 3:
        return "loss"
    return "neutral"


def parse_reference(path, cn_col, neutral_cn):
    """SNP-array BED -> list of events {chrom,start,end,length,cn,dir,label}."""
    events = []
    skipped = 0
    idx = cn_col - 1                       # 1-based -> 0-based
    with open(path) as fh:
        for line in fh:
            raw = line.rstrip("\n")
            if not raw.strip() or raw.startswith(("#", "track", "browser")):
                continue
            f = raw.split("\t")
            if len(f) <= idx:
                continue
            try:
                start, end = int(f[1]), int(f[2])
                cn = float(f[idx])
            except ValueError:
                skipped += 1               # header / malformed -> skip
                continue
            d = cn_to_dir(cn, neutral_cn)
            if d == "neutral":
                continue
            label = f[3] if len(f) > 3 else ""     # e.g. loss / amplification
            events.append({"chrom": norm_chrom(f[0]), "start": start, "end": end,
                           "length": end - start, "cn": cn, "dir": d, "label": label})
    if skipped:
        print(f"[overlap] reference: skipped {skipped} header/malformed line(s)")
    return events


def parse_query(path, group=None):
    """inferCNV pred_cnv_regions.dat -> {dir: {chrom: [(start,end), ...]}}.
    Detects the header if present; otherwise assumes the fixed column order
    cell_group_name, cnv_name, state, chr, start, end."""
    by_dir = {"gain": defaultdict(list), "loss": defaultdict(list)}
    groups = defaultdict(int)
    order = {"cell_group_name": 0, "state": 2, "chr": 3, "start": 4, "end": 5}

    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not lines:
        sys.exit(f"[overlap] query file is empty: {path}")

    first = lines[0].split("\t")
    has_header = "state" in first and "chr" in first and "start" in first
    if has_header:
        col = {name: first.index(name) for name in
               ("cell_group_name", "state", "chr", "start", "end") if name in first}
        for name in ("state", "chr", "start", "end"):
            if name not in col:
                sys.exit(f"[overlap] query header missing required column '{name}': {first}")
        data = lines[1:]
    else:
        col = order
        data = lines

    n_kept = n_other_group = 0
    for ln in data:
        f = ln.split("\t")
        if len(f) <= col["end"]:
            continue
        grp = f[col["cell_group_name"]] if "cell_group_name" in col and len(f) > col["cell_group_name"] else "NA"
        groups[grp] += 1
        if group is not None and grp != group:
            n_other_group += 1
            continue
        try:
            state = int(float(f[col["state"]]))
            start, end = int(f[col["start"]]), int(f[col["end"]])
        except ValueError:
            continue
        d = state_to_dir(state)
        if d == "neutral":
            continue
        by_dir[d][norm_chrom(f[col["chr"]])].append((start, end))
        n_kept += 1

    print(f"[overlap] query cell groups found: " +
          ", ".join(f"{g} ({c})" for g, c in sorted(groups.items())))
    if group is not None:
        print(f"[overlap] restricted to group '{group}': kept {n_kept} region(s), "
              f"ignored {n_other_group} from other groups")
    else:
        print(f"[overlap] using ALL groups (union): {n_kept} non-neutral region(s)")

    for d in by_dir:
        for chrom in by_dir[d]:
            by_dir[d][chrom] = merge_intervals(by_dir[d][chrom])
    return by_dir


# --------------------------------------------------------------------------- #
# interval math
# --------------------------------------------------------------------------- #
def merge_intervals(intervals):
    """Merge overlapping/touching [start,end) intervals; returns sorted list."""
    if not intervals:
        return []
    ivs = sorted(intervals)
    merged = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(m) for m in merged]


def overlap_bp(start, end, merged):
    """Total bp of [start,end) covered by a sorted, merged interval list."""
    total = 0
    for s, e in merged:
        if e <= start:
            continue
        if s >= end:
            break
        total += min(end, e) - max(start, s)
    return total


# --------------------------------------------------------------------------- #
# comparison
# --------------------------------------------------------------------------- #
def compare(events, query, min_overlap):
    opp = {"gain": "loss", "loss": "gain"}
    rows = []
    for ev in events:
        same = overlap_bp(ev["start"], ev["end"],
                          query[ev["dir"]].get(ev["chrom"], []))
        contra = overlap_bp(ev["start"], ev["end"],
                            query[opp[ev["dir"]]].get(ev["chrom"], []))
        frac_same = same / ev["length"] if ev["length"] else 0.0
        frac_contra = contra / ev["length"] if ev["length"] else 0.0
        rows.append({**ev, "same_bp": same, "same_frac": frac_same,
                     "contra_bp": contra, "contra_frac": frac_contra,
                     "matched": frac_same >= min_overlap})
    return rows


def summarize(rows, min_overlap):
    """Aggregate overall and per direction."""
    out = {}
    for key in ("all", "gain", "loss"):
        sel = rows if key == "all" else [r for r in rows if r["dir"] == key]
        n = len(sel)
        n_match = sum(r["matched"] for r in sel)
        ref_bp = sum(r["length"] for r in sel)
        same_bp = sum(r["same_bp"] for r in sel)
        contra_bp = sum(r["contra_bp"] for r in sel)
        out[key] = {
            "n_events": n,
            "n_matched": n_match,
            "pct_events_matched": 100.0 * n_match / n if n else float("nan"),
            "ref_bp": ref_bp,
            "same_bp": same_bp,
            "pct_bp_matched": 100.0 * same_bp / ref_bp if ref_bp else float("nan"),
            "contra_bp": contra_bp,
            "pct_bp_contradicted": 100.0 * contra_bp / ref_bp if ref_bp else float("nan"),
        }
    return out


# --------------------------------------------------------------------------- #
# reporting
# --------------------------------------------------------------------------- #
def print_report(summ, rows, min_overlap):
    def pct(x):
        return "  nan" if x != x else f"{x:5.1f}"

    print()
    print("=" * 72)
    print(f"SNP array (reference) vs inferCNV (query)   [match >= {min_overlap:.0%} "
          "same-direction overlap of each reference event]")
    print("=" * 72)
    hdr = f"{'direction':<8} {'events':>7} {'matched':>8} {'%evt':>6}   " \
          f"{'ref_bp':>13} {'%bp_match':>10} {'%bp_contra':>11}"
    print(hdr)
    print("-" * len(hdr))
    for key in ("all", "gain", "loss"):
        s = summ[key]
        print(f"{key:<8} {s['n_events']:>7} {s['n_matched']:>8} {pct(s['pct_events_matched'])}   "
              f"{s['ref_bp']:>13,} {pct(s['pct_bp_matched'])}      {pct(s['pct_bp_contradicted'])}")
    print()
    print("  % evt        = Metric 1: reference events matched / all reference events")
    print("  % bp_match   = Metric 2: same-direction overlapping bp / all reference-event bp")
    print("                 (denominator includes non-matching events)")
    print("  % bp_contra  = reference bp overlapping the OPPOSITE inferCNV call (disagreement)")

    # sensitivity of the event-match metric to the overlap threshold
    print()
    print("  event-match % vs overlap threshold (overall):")
    line = "   "
    for thr in (0.10, 0.25, 0.50, 0.75, 0.90):
        n = len(rows)
        m = sum(r["same_frac"] >= thr for r in rows)
        line += f"  >={thr:.0%}: {100.0*m/n if n else float('nan'):4.0f}%"
    print(line)
    print()


def write_tables(rows, summ, prefix, min_overlap):
    per_event = prefix + ".per_event.tsv"
    with open(per_event, "w") as fh:
        fh.write("chrom\tstart\tend\tlength\tcn\tlabel\tdirection\t"
                 "same_overlap_bp\tsame_overlap_frac\t"
                 "contra_overlap_bp\tcontra_overlap_frac\tmatched\n")
        for r in sorted(rows, key=lambda r: (r["chrom"], r["start"])):
            fh.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['length']}\t"
                     f"{r['cn']}\t{r['label']}\t{r['dir']}\t"
                     f"{r['same_bp']}\t{r['same_frac']:.4f}\t"
                     f"{r['contra_bp']}\t{r['contra_frac']:.4f}\t"
                     f"{int(r['matched'])}\n")

    summary = prefix + ".summary.tsv"
    with open(summary, "w") as fh:
        fh.write("direction\tn_events\tn_matched\tpct_events_matched\t"
                 "ref_bp\tsame_bp\tpct_bp_matched\tcontra_bp\tpct_bp_contradicted\n")
        for key in ("all", "gain", "loss"):
            s = summ[key]
            fh.write(f"{key}\t{s['n_events']}\t{s['n_matched']}\t"
                     f"{s['pct_events_matched']:.2f}\t{s['ref_bp']}\t{s['same_bp']}\t"
                     f"{s['pct_bp_matched']:.2f}\t{s['contra_bp']}\t"
                     f"{s['pct_bp_contradicted']:.2f}\n")
    print(f"[overlap] wrote {per_event}")
    print(f"[overlap] wrote {summary}")


# --------------------------------------------------------------------------- #
def run(reference, query_path, cn_col=13, min_overlap=0.5, group=None,
        neutral_cn=2.0, out_prefix=None):
    events = parse_reference(reference, cn_col, neutral_cn)
    if not events:
        sys.exit("[overlap] no gain/loss events in the reference — nothing to compare.")
    n_gain = sum(e["dir"] == "gain" for e in events)
    n_loss = sum(e["dir"] == "loss" for e in events)
    print(f"[overlap] reference events: {len(events)} ({n_gain} gain, {n_loss} loss)")

    query = parse_query(query_path, group)
    rows = compare(events, query, min_overlap)
    summ = summarize(rows, min_overlap)
    print_report(summ, rows, min_overlap)
    if out_prefix:
        write_tables(rows, summ, out_prefix, min_overlap)
    return summ, rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-r", "--reference", required=True,
                   help="SNP-array CNV BED on hg38 (lift with liftover_snp_array.py first)")
    p.add_argument("-q", "--query", required=True,
                   help="inferCNV pred_cnv_regions.dat")
    p.add_argument("--cn-col", type=int, default=13,
                   help="1-based copy-number column in the reference BED (default 13)")
    p.add_argument("--min-overlap", type=float, default=0.5,
                   help="min same-direction overlap fraction of a reference event to "
                        "count it as matched (default 0.5)")
    p.add_argument("--neutral-cn", type=float, default=2.0,
                   help="copy number treated as diploid/neutral (default 2.0)")
    p.add_argument("--group", default=None,
                   help="restrict inferCNV to one cell_group_name (default: union of all)")
    p.add_argument("--out-prefix", default=None,
                   help="write <prefix>.per_event.tsv and <prefix>.summary.tsv")
    args = p.parse_args(argv)
    run(args.reference, args.query, args.cn_col, args.min_overlap, args.group,
        args.neutral_cn, args.out_prefix)


if __name__ == "__main__":
    main()
