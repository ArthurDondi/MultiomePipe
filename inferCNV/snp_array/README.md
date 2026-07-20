# Compare inferCNV against cell-line SNP-array CNVs

Validate the inferCNV copy-number calls (from `inferCNV/run_infercnv.R`) against
SNP-array CNV segments from the matching cell lines, taking the **SNP array as
the reference / ground truth**.

```
SNP-array CNV BED (hg19)                 inferCNV results/<sample>/
        |                                17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat
        | 1) liftover_snp_array.py                    |
        v                                             |
  cellline.hg38.bed  --------- 2) compare_cnv_overlap.py --------->  % events & % bp matching
```

Two steps because the two files are on **different genome builds**: the SNP array
is GRCh37/hg19 (e.g. a chr6 loss ending at 171,051,066 sits past the end of hg38
chr6 (170,805,979) but inside hg19 chr6 (171,115,067); the MYCN amp at
chr2:15.0-16.4 Mb matches hg19 MYCN), while inferCNV runs on the hg38 gencode v50
gene ordering built by `export_for_infercnv.py`. So the array is lifted to hg38
first, then overlapped.

## Install

Only the liftover step has an extra dependency:

```bash
pip install pyliftover          # or: conda install -c bioconda pyliftover
```

`compare_cnv_overlap.py` is pure standard-library Python (no pandas/numpy needed).

## 1. Lift the SNP array hg19 -> hg38

```bash
python inferCNV/snp_array/liftover_snp_array.py \
    -i cellline_SKNBE2c.hg19.bed \
    -o cellline_SKNBE2c.hg38.bed
```

On first use pyliftover downloads UCSC's `hg19ToHg38.over.chain.gz`. On a cluster
with no outbound internet, fetch the chain once and pass it explicitly:

```bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
python inferCNV/snp_array/liftover_snp_array.py -i in.hg19.bed -o out.hg38.bed \
    --chain hg19ToHg38.over.chain.gz
```

Every original column is preserved; only the coordinate columns are rewritten
(cols 1-3, plus the duplicate start/end/length cols 7/8/11 when they mirror the
originals). Segments whose endpoints don't map, land on two different hg38
chromosomes, or change length by more than `--max-len-change` (default 25%) are
dropped to `<output>.unmapped` with a reason column.

## 2. Overlap against inferCNV

```bash
python inferCNV/snp_array/compare_cnv_overlap.py \
    -r cellline_SKNBE2c.hg38.bed \
    -q results/BMO-SKNBE2c/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat \
    --out-prefix results/BMO-SKNBE2c/snp_vs_infercnv
```

Point `-q` at the HMM regions file inferCNV writes when `USE_HMM = TRUE`. Prefer
the Bayes-filtered `HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.5.pred_cnv_regions.dat`
if present (cleaner); the `17_HMM_pred...` file is the pre-filter version and has
the identical layout. Match each SNP-array cell line to its `results/<sample>/`
folder (e.g. SK-N-BE(2)c -> `BMO-SKNBE2c`).

### What it computes

Only the **direction** of each event is compared (never absolute copy number):

| source | column | rule |
|--------|--------|------|
| SNP array (reference) | copy number, default **col 13** (`--cn-col`) | `>2` gain, `<2` loss, `==2` neutral (ignored) |
| inferCNV (query) | HMM `state` (1-6, 3 = diploid) | `>3` gain, `<3` loss, `==3` neutral (never emitted) |

A reference event **matches** when same-direction inferCNV calls cover at least
`--min-overlap` (default **0.50**) of its length. Overlapping inferCNV intervals
are merged per direction first, so base pairs are never double-counted (this also
unions several `cell_group_name` subclusters cleanly; use `--group` to restrict to
one).

Two headline numbers (SNP array as reference), reported overall and split into
gains / losses:

1. **% events matching** = matched reference events / all reference events.
2. **% bp matching** = same-direction overlapping bp / all reference-event bp.
   The denominator includes reference events that never clear the per-event
   threshold, so this is a pure base-pair directional concordance.

Opposite-direction (`% bp contradicted`) bp are also reported as QC, along with a
small table of how the event-match rate moves as the overlap threshold changes.

Outputs: `<prefix>.summary.tsv` (the table above) and `<prefix>.per_event.tsv`
(every reference event with its same/opposite overlap and matched flag, for
inspecting misses).

## Caveats worth remembering

- **Relative vs absolute / reference-less runs.** inferCNV reports CN *relative to
  its baseline*. Pure cell-line samples with no normal cells run reference-less
  (baseline = the cell line's own mean), so uniform whole-chromosome gains or a
  genome doubling can read as "neutral" and lower the apparent concordance. For
  the fairest comparison, run those samples against a real diploid reference
  (borrow normal T/NK/B cells, or use `MODE = "joint"`).
- **Direction only.** A CN=40 amplification and a CN=3 gain are both just "gain"
  here; this measures agreement on gain/loss direction, not magnitude.
- **Resolution.** inferCNV works off gene windows, so focal / sub-Mb array calls
  (e.g. the intragenic IGF2BP3 deletion, the MIR129-1 gain) are usually invisible
  and will show as unmatched — that is expected, not a bug. Large arm-level events
  (2p/MYCN, 17q, 1p, 6q) are where agreement should appear.
- **Build.** Both files must be hg38 before step 2; don't skip the liftover.
