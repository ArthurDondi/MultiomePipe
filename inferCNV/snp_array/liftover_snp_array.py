#!/usr/bin/env python
"""
Lift a cell-line SNP-array CNV BED from hg19/GRCh37 to hg38/GRCh38 so it lines up
with inferCNV output (inferCNV runs on the hg38 gencode v50 gene ordering built
by export_for_infercnv.py from the CellRanger reference).

Why this is needed: the SNP-array segments are on GRCh37. For example a chr6 loss
ending at 171,051,066 sits *past* the end of hg38 chr6 (170,805,979) but inside
hg19 chr6 (171,115,067); a MYCN amplification at chr2:15.0-16.4 Mb also matches
hg19 MYCN (~16.08 Mb). Overlapping the two builds directly would be meaningless,
so the array has to be lifted to hg38 first.

Engine: pyliftover (pure Python) ->  pip install pyliftover
On first use pyliftover downloads UCSC's hg19ToHg38.over.chain.gz. On an offline
cluster, fetch the chain once and pass it with --chain:
  http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

Each segment's start and end are lifted independently. If an endpoint lands in a
chain gap (common in sub-telomeric / peri-centromeric sequence) it is walked
INWARD up to --max-nudge bp until a base maps, so a whole-arm CNV is not dropped
just because its boundary base has no hg38 counterpart (a few kb off a Mb-scale
segment is negligible; the trim is reported). A segment is only *dropped* (written
to <output>.unmapped instead) when
  - an endpoint still fails to map after nudging up to --max-nudge bp,
  - the two endpoints land on different hg38 chromosomes (segment straddles a
    build rearrangement), or
  - the lifted length differs from the original by more than --max-len-change.
Every original column is preserved; only the coordinate columns are rewritten
(cols 1-3, and the duplicate start/end/length columns 7/8/11 when they mirror the
originals, as they do in the sample file), so the output is a drop-in hg38 BED for
compare_cnv_overlap.py.

Usage:
    python liftover_snp_array.py -i cellline.hg19.bed -o cellline.hg38.bed
    python liftover_snp_array.py -i in.bed -o out.bed --chain hg19ToHg38.over.chain.gz
"""

import os
import sys
import argparse


def load_liftover(chain, from_build, to_build):
    try:
        from pyliftover import LiftOver
    except ImportError:
        sys.exit("[liftover] pyliftover not installed. Run:  pip install pyliftover\n"
                 "          (or conda install -c bioconda pyliftover)")
    if chain:
        if not os.path.exists(chain):
            sys.exit(f"[liftover] --chain file not found: {chain}")
        print(f"[liftover] using local chain file: {chain}")
        return LiftOver(chain)

    # Auto-download path: fine on a networked machine, but on a locked-down
    # cluster the fetch silently yields no file and pyliftover then blows up deep
    # inside with an unhelpful AttributeError. Catch that and point at --chain.
    chain_name = f"{from_build}To{to_build[0].upper()}{to_build[1:]}.over.chain.gz"
    chain_url = (f"http://hgdownload.soe.ucsc.edu/goldenPath/{from_build}/"
                 f"liftOver/{chain_name}")
    print(f"[liftover] fetching the UCSC {from_build}->{to_build} chain via pyliftover "
          "(pass --chain to work offline)")
    try:
        lo = LiftOver(from_build, to_build)
    except Exception as e:
        lo = None
        err = f"{type(e).__name__}: {e}"
    else:
        err = None
    if lo is None or getattr(lo, "chain_file", None) is None:
        sys.exit(
            f"[liftover] could not download the {from_build}->{to_build} chain "
            f"automatically{(' (' + err + ')') if err else ''}.\n"
            f"[liftover] Download it once, then re-run with --chain:\n"
            f"[liftover]     wget {chain_url}\n"
            f"[liftover]     python liftover_snp_array.py -i IN -o OUT --chain {chain_name}")
    return lo


def norm_chrom(c):
    """UCSC chain files are chr-prefixed; make sure the query is too."""
    return c if c.startswith("chr") else "chr" + c


def lift_point(lo, chrom, pos):
    """Lift one 0-based coordinate. Returns (chrom, pos) or None."""
    res = lo.convert_coordinate(chrom, pos)
    if not res:                       # None (chrom absent) or [] (unmapped)
        return None
    tgt_chrom, tgt_pos = res[0][0], res[0][1]
    return tgt_chrom, tgt_pos


def lift_point_rescue(lo, chrom, pos, inward, max_nudge, step=1000):
    """Lift `pos`; if it falls in a chain gap, walk INWARD until a base maps.
    `inward` is +1 for a start (move right, into the segment) or -1 for an end
    (move left). Steps by `step` bp up to `max_nudge`. Returns
    (chrom, mapped_pos, shift_bp) — shift_bp is how far we had to move (0 if the
    original base mapped) — or None if nothing maps within max_nudge."""
    r = lift_point(lo, chrom, pos)
    if r is not None:
        return r[0], r[1], 0
    moved = step
    while moved <= max_nudge:
        r = lift_point(lo, chrom, pos + inward * moved)
        if r is not None:
            return r[0], r[1], moved
        moved += step
    return None


def lift_interval(lo, chrom, start, end, max_nudge=0, step=1000):
    """Lift a BED interval [start, end). end is exclusive, so the last included
    base is end-1; we lift that and add 1 back. When max_nudge > 0, an unmappable
    endpoint is walked inward (start right, end left) up to min(max_nudge, ~half
    the segment) bp so the segment is not dropped for a single unmappable boundary
    base. Returns (chrom, new_start, new_end, trim_bp) on success — trim_bp is the
    hg19 bp shaved off the ends by rescue (0 if none) — or (None, reason)."""
    # cap the nudge so start and end can never cross past each other
    cap = min(max_nudge, max(0, (end - start) // 2 - 1)) if max_nudge else 0
    a = lift_point_rescue(lo, chrom, start, +1, cap, step)      # start moves right
    b = lift_point_rescue(lo, chrom, end - 1, -1, cap, step)    # end moves left
    if a is None or b is None:
        return None, "unmapped_endpoint"
    if a[0] != b[0]:
        return None, "endpoints_on_different_chromosomes"
    lo_pos, hi_pos = sorted((a[1], b[1]))   # handles +/- strand mapping
    return a[0], lo_pos, hi_pos + 1, a[2] + b[2]


def run(input_bed, output_bed, chain=None, from_build="hg19", to_build="hg38",
        max_len_change=0.25, max_nudge=200000):
    lo = load_liftover(chain, from_build, to_build)

    unmapped_path = output_bed + ".unmapped"
    n_in = n_out = n_skip_hdr = 0
    reasons = {}
    rescued = []                       # (chrom, start, end, trim_bp)

    with open(input_bed) as fin, \
         open(output_bed, "w") as fout, \
         open(unmapped_path, "w") as funmap:
        for line in fin:
            raw = line.rstrip("\n")
            if not raw.strip() or raw.startswith(("#", "track", "browser")):
                continue
            f = raw.split("\t")
            if len(f) < 3:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                n_skip_hdr += 1          # header row (non-integer coords) -> skip
                continue
            n_in += 1
            chrom = norm_chrom(f[0])
            old_len = end - start

            res = lift_interval(lo, chrom, start, end, max_nudge)
            if res[0] is None:
                reasons[res[1]] = reasons.get(res[1], 0) + 1
                funmap.write(raw + "\t" + res[1] + "\n")
                continue

            new_chrom, new_start, new_end, trim = res
            new_len = new_end - new_start
            if old_len > 0 and abs(new_len - old_len) / old_len > max_len_change:
                reasons["length_change"] = reasons.get("length_change", 0) + 1
                funmap.write(raw + f"\tlength_change({old_len}->{new_len})\n")
                continue
            if trim > 0:
                rescued.append((chrom, start, end, trim))

            # rewrite coordinate columns; also fix the duplicate start/end (cols
            # 7/8) and length (col 11) *only* when they mirror the hg19 values.
            f[0], f[1], f[2] = new_chrom, str(new_start), str(new_end)
            if len(f) >= 8 and f[6] == str(start) and f[7] == str(end):
                f[6], f[7] = str(new_start), str(new_end)
            if len(f) >= 11 and f[10] == str(old_len):
                f[10] = str(new_len)
            fout.write("\t".join(f) + "\n")
            n_out += 1

    print(f"[liftover] input events        : {n_in}"
          + (f"  (+{n_skip_hdr} header/comment line(s) skipped)" if n_skip_hdr else ""))
    print(f"[liftover] lifted to {to_build}      : {n_out}")
    if rescued:
        print(f"[liftover] rescued (nudged)    : {len(rescued)}  "
              f"(unmappable boundary, walked inward <= {max_nudge:,} bp)")
        for c, s, e, tr in rescued:
            print(f"[liftover]     - {c}:{s:,}-{e:,} trimmed {tr:,} bp")
    n_unmapped = n_in - n_out
    print(f"[liftover] dropped (unmapped)  : {n_unmapped}")
    for r, c in sorted(reasons.items()):
        print(f"[liftover]     - {r}: {c}")
    print(f"[liftover] wrote {output_bed}")
    if n_unmapped:
        print(f"[liftover] dropped rows -> {unmapped_path} (with a reason column)")
    else:
        os.remove(unmapped_path)
    return output_bed


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="SNP-array CNV BED on hg19")
    p.add_argument("-o", "--output", required=True, help="output BED on hg38")
    p.add_argument("--chain", default=None,
                   help="local UCSC over.chain(.gz) (default: pyliftover downloads it)")
    p.add_argument("--from-build", default="hg19")
    p.add_argument("--to-build", default="hg38")
    p.add_argument("--max-len-change", type=float, default=0.25,
                   help="drop a segment if |lifted-len - orig-len|/orig-len exceeds "
                        "this (default 0.25 = 25%%; guards against build rearrangements)")
    p.add_argument("--max-nudge", type=int, default=200000,
                   help="if an endpoint is unmappable, walk it inward up to this many "
                        "bp to rescue the segment (default 200000; 0 disables). A few "
                        "kb off a Mb-scale CNV is negligible; the trim is reported.")
    args = p.parse_args(argv)
    run(args.input, args.output, args.chain, args.from_build, args.to_build,
        args.max_len_change, args.max_nudge)


if __name__ == "__main__":
    main()