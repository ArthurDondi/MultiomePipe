#!/usr/bin/env python
"""
Overlap a cell-line SNP-array CNV BED (hg38) against an inferCNV HMM prediction
file, taking the SNP array as the reference / ground truth.

    reference = SNP array  (pred_cnv_regions-style truth for the cell line)
    query     = inferCNV   (17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat)

Both files must be on the SAME genome build. The array is hg19, so lift it first
with liftover_snp_array.py; inferCNV is already hg38.

Direction of an event (only the direction is compared, never absolute copies):
  reference (SNP array): category column, default col 4 (--type-col)
        gain / wc_gain / amplification -> gain ;  loss / wc_loss -> loss ;
        anything else -> the segment is ignored  (wc = whole chromosome)
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

import os
import sys
import argparse
from collections import defaultdict


# --------------------------------------------------------------------------- #
# parsing
# --------------------------------------------------------------------------- #
def norm_chrom(c):
    return c if c.startswith("chr") else "chr" + c


def label_to_dir(label):
    """Direction from the SNP-array category column (default col 4): any label
    containing 'gain' or 'amplification' (gain, wc_gain, small_gain, amplification;
    wc_ = whole chromosome) -> gain; containing 'loss' (loss, wc_loss, small_loss)
    -> loss; anything else -> None so the segment is ignored. Case-insensitive
    substring match."""
    s = label.lower()
    if "loss" in s:
        return "loss"
    if "gain" in s or "amplification" in s:
        return "gain"
    return None


def state_to_dir(state, neutral_state=3):
    """inferCNV HMM state -> direction, relative to the model's diploid state.
    i6 model: states 1,2 loss / 3 neutral / 4,5,6 gain  (neutral_state=3).
    i3 model: state  1   loss / 2 neutral / 3     gain  (neutral_state=2)."""
    if state > neutral_state:
        return "gain"
    if state < neutral_state:
        return "loss"
    return "neutral"


def neutral_state_for(query_path, hmm_i=None):
    """Return (neutral_state, hmm_i). If hmm_i is None, autodetect the HMM model
    from the pred_cnv_regions filename ('...predHMMi3...' or '...HMMi6...'); default
    to i6 (inferCNV's own default) if nothing matches. Neutral is the middle state:
    i3 -> 2, i6 -> 3."""
    if hmm_i is None:
        base = os.path.basename(query_path)
        if "HMMi3" in base:
            hmm_i = 3
        elif "HMMi6" in base:
            hmm_i = 6
        else:
            hmm_i = 6
    if hmm_i not in (3, 6):
        sys.exit(f"[overlap] --hmm-i must be 3 or 6, got {hmm_i}")
    return {3: 2, 6: 3}[hmm_i], hmm_i


def parse_reference(path, type_col=4, cn_col=13):
    """SNP-array BED -> list of events {chrom,start,end,length,cn,dir,label}.

    Direction comes from the category column `type_col` (default col 4) via
    label_to_dir: gain/wc_gain/amplification -> gain, loss/wc_loss -> loss; a
    segment whose category has none of those terms is IGNORED. `cn_col` (default 13)
    is read best-effort only for the recorded copy number (blank if non-numeric or
    absent); it no longer drives the direction, so a wrong/empty cn column can't
    zero out the reference."""
    events = []
    ignored = 0                            # no gain/loss/amplification category
    malformed = 0                          # bad coordinates
    ti = type_col - 1                      # 1-based -> 0-based
    ci = cn_col - 1
    with open(path) as fh:
        for line in fh:
            raw = line.rstrip("\n")
            if not raw.strip() or raw.startswith(("#", "track", "browser")):
                continue
            f = raw.split("\t")
            if len(f) < 3:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                malformed += 1             # header / non-integer coordinates
                continue
            label = f[ti] if len(f) > ti else ""
            d = label_to_dir(label)
            if d is None:                  # no gain/loss category in col type_col
                ignored += 1
                continue
            cn = None
            if 0 <= ci < len(f):
                try:
                    cn = float(f[ci])
                except ValueError:
                    cn = None
            events.append({"chrom": norm_chrom(f[0]), "start": start, "end": end,
                           "length": end - start, "cn": cn, "dir": d, "label": label})
    if malformed:
        print(f"[overlap] reference: skipped {malformed} line(s) with non-integer coords")
    if ignored:
        print(f"[overlap] reference: ignored {ignored} line(s) with no gain/loss/"
              f"amplification category in col {type_col}")
    return events


def parse_query_by_group(path, group=None, quiet=False, neutral_state=3):
    """inferCNV pred_cnv_regions.dat -> {cell_group_name: {dir: {chrom: [merged]}}}.

    Each cell_group_name is one inferCNV group: a cell type in analysis_mode
    "samples", or a subclone/clone in "subclusters" mode (e.g. 'neuron.neuron_s8').
    `neutral_state` is the model's diploid state (i6 -> 3, i3 -> 2; see
    neutral_state_for). Detects the header if present; otherwise assumes the fixed
    column order cell_group_name, cnv_name, state, chr, start, end. If `group` is
    given, only that clone is kept."""
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

    raw = defaultdict(lambda: {"gain": defaultdict(list), "loss": defaultdict(list)})
    counts = defaultdict(int)
    for ln in data:
        f = ln.split("\t")
        if len(f) <= col["end"]:
            continue
        gi = col.get("cell_group_name", 0)
        grp = f[gi] if len(f) > gi else "NA"
        counts[grp] += 1
        if group is not None and grp != group:
            continue
        try:
            state = int(float(f[col["state"]]))
            start, end = int(f[col["start"]]), int(f[col["end"]])
        except ValueError:
            continue
        d = state_to_dir(state, neutral_state)
        if d == "neutral":
            continue
        raw[grp][d][norm_chrom(f[col["chr"]])].append((start, end))

    for grp in raw:
        for d in ("gain", "loss"):
            for chrom in raw[grp][d]:
                raw[grp][d][chrom] = merge_intervals(raw[grp][d][chrom])

    if not quiet:
        print("[overlap] query groups found: " +
              ", ".join(f"{g} ({c})" for g, c in sorted(counts.items())))
    if group is not None and group not in counts:
        sys.exit(f"[overlap] --group '{group}' not found. Available: {sorted(counts)}")
    return dict(raw)


def union_query(by_group):
    """Collapse {group: by_dir} into one by_dir unioned (merged) across clones."""
    union = {"gain": defaultdict(list), "loss": defaultdict(list)}
    for bd in by_group.values():
        for d in ("gain", "loss"):
            for chrom, ivs in bd[d].items():
                union[d][chrom].extend(ivs)
    for d in ("gain", "loss"):
        for chrom in union[d]:
            union[d][chrom] = merge_intervals(union[d][chrom])
    return union


def parse_query(path, group=None, neutral_state=3):
    """inferCNV pred_cnv_regions.dat -> single {dir: {chrom: [merged]}}, unioned
    across all clones (or restricted to one with `group`)."""
    by_group = parse_query_by_group(path, group, neutral_state=neutral_state)
    if group is not None:
        print(f"[overlap] restricted to group '{group}'")
    else:
        print(f"[overlap] using ALL groups (union of {len(by_group)} group(s))")
    return union_query(by_group)


def parse_groupings(path, quiet=False):
    """inferCNV observation_groupings.txt -> {subcluster_name: n_cells}.

    R write.table layout: a header row of quoted column names, then one row per
    cell: "<barcode>" "<Dendrogram Group>" "<color>" "<Annotation Group>" "<color>".
    The 2nd quoted field (Dendrogram Group) is the subcluster id (e.g. neuron_s8),
    which matches the suffix of a pred_cnv_regions cell_group_name."""
    import re
    counts = defaultdict(int)
    n = 0
    with open(path) as fh:
        for ln in fh:
            toks = re.findall(r'"([^"]*)"', ln) or ln.split()
            if not toks or "Dendrogram Group" in toks or len(toks) < 2:
                continue                       # header / blank / malformed
            counts[toks[1]] += 1
            n += 1
    if not counts:
        print(f"[overlap] WARNING: no clone sizes parsed from {path}")
    elif not quiet:
        print(f"[overlap] groupings: {n} cells across {len(counts)} subcluster(s)")
    return dict(counts)


def clone_size(cell_group_name, sizes):
    """Resolve a pred_cnv_regions cell_group_name (e.g. 'neuron.neuron_s8') to a
    cell count from parse_groupings (keyed by subcluster, e.g. 'neuron_s8'): try
    exact, then the part after the first '.', then the longest subcluster key the
    name ends with. Returns None if unresolved."""
    if not sizes:
        return None
    if cell_group_name in sizes:
        return sizes[cell_group_name]
    tail = cell_group_name.split(".", 1)[1] if "." in cell_group_name else cell_group_name
    if tail in sizes:
        return sizes[tail]
    cands = [k for k in sizes if cell_group_name.endswith(k)]
    return sizes[max(cands, key=len)] if cands else None


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
            cn = "" if r["cn"] is None else r["cn"]
            fh.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['length']}\t"
                     f"{cn}\t{r['label']}\t{r['dir']}\t"
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
def run(reference, query_path, type_col=4, cn_col=13, min_overlap=0.5, group=None,
        out_prefix=None, hmm_i=None):
    events = parse_reference(reference, type_col, cn_col)
    if not events:
        sys.exit("[overlap] no gain/loss events in the reference — nothing to compare.")
    n_gain = sum(e["dir"] == "gain" for e in events)
    n_loss = sum(e["dir"] == "loss" for e in events)
    print(f"[overlap] reference events: {len(events)} ({n_gain} gain, {n_loss} loss)")

    neutral_state, hmm_i = neutral_state_for(query_path, hmm_i)
    print(f"[overlap] inferCNV model: HMMi{hmm_i} (neutral = state {neutral_state})")
    query = parse_query(query_path, group, neutral_state)
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
    p.add_argument("--type-col", type=int, default=4,
                   help="1-based category column giving the event direction "
                        "(gain/wc_gain/amplification vs loss/wc_loss); default 4")
    p.add_argument("--cn-col", type=int, default=13,
                   help="1-based copy-number column, recorded for reference only "
                        "(default 13); blank if non-numeric — does not affect direction")
    p.add_argument("--min-overlap", type=float, default=0.5,
                   help="min same-direction overlap fraction of a reference event to "
                        "count it as matched (default 0.5)")
    p.add_argument("--group", default=None,
                   help="restrict inferCNV to one cell_group_name (default: union of all)")
    p.add_argument("--hmm-i", type=int, choices=(3, 6), default=None,
                   help="inferCNV HMM model: i3 (1 loss/2 neutral/3 gain) or i6 "
                        "(1,2 loss/3 neutral/4-6 gain). Default: autodetect from the "
                        "query filename (predHMMi3 / HMMi6), else i6.")
    p.add_argument("--out-prefix", default=None,
                   help="write <prefix>.per_event.tsv and <prefix>.summary.tsv")
    args = p.parse_args(argv)
    run(args.reference, args.query, args.type_col, args.cn_col, args.min_overlap,
        args.group, args.out_prefix, args.hmm_i)


if __name__ == "__main__":
    main()
