# Compare inferCNV against cell-line SNP-array CNVs

Validate the inferCNV copy-number calls (from `inferCNV/run_infercnv.R`) against
SNP-array CNV segments from the matching cell lines, taking the **SNP array as
the reference / ground truth**.

```
SNP-array CNV BED (hg19)                 inferCNV results/<sample>/
        |                                *.hmm_mode-{samples,subclusters}.pred_cnv_regions.dat
        | 1) liftover_snp_array.py                    |  (+ observation_groupings.txt for clone sizes)
        v                                             |
  cellline.hg38.bed  --------- 2) compare_cnv_overlap.py --------->  % events & % bp matching
        |                                             |
        +----------- 3) plot_cnv_overlay.py ----------+--------->  genome-wide overlay (PNG/PDF)
        |                                             |
        +----------- 4) plot_clones.py ---------------+--------->  per-clone overlays + clone summary
                     (analysis_mode = "subclusters")             (each clone = one query)
```

Two steps because the two files are on **different genome builds**: the SNP array
is GRCh37/hg19 (e.g. a chr6 loss ending at 171,051,066 sits past the end of hg38
chr6 (170,805,979) but inside hg19 chr6 (171,115,067); the MYCN amp at
chr2:15.0-16.4 Mb matches hg19 MYCN), while inferCNV runs on the hg38 gencode v50
gene ordering built by `export_for_infercnv.py`. So the array is lifted to hg38
first, then overlapped.

## Install

```bash
pip install pyliftover          # step 1 only (or: conda install -c bioconda pyliftover)
pip install matplotlib          # step 3 only
```

`compare_cnv_overlap.py` (step 2) is pure standard-library Python (no deps). The
`scverse` conda env already has matplotlib.

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
originals). Segments whose endpoints can't be placed on the source chromosome, or
that change length by more than `--max-len-change` (default 25%), are dropped to
`<output>.unmapped` with a reason column.

**Endpoint rescue.** hg19->hg38 keeps a segment on the same chromosome, so both
endpoints are required to land on the source chromosome. A boundary base in a
chain gap, or one that maps to a *different* chromosome (sub-telomeric /
peri-centromeric repeats do both — e.g. a whole-arm chr1 gain whose telomeric end
maps onto another chromosome, or a 91 Mb chr12 gain starting at the chr12:149,960
sub-telomere), would otherwise sink a whole-arm CNV. So a bad endpoint is walked
**inward** in 1 kb steps up to `--max-nudge` bp (default 200 kb; `0` disables)
until a base maps to the source chromosome; a few kb trimmed off a Mb-scale segment
is negligible and every rescued segment is reported with its trim. The nudge is
capped at half the segment so a tiny segment with a deep-gap boundary still drops
rather than collapsing.

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
| SNP array (reference) | **category**, default **col 4** (`--type-col`) | `gain` / `wc_gain` / `amplification` -> gain; `loss` / `wc_loss` -> loss; else ignored |
| inferCNV (query) | HMM `state` | `> neutral` gain, `< neutral` loss, `== neutral` (never emitted) |

**Reference direction (`--type-col`, default col 4).** The array caller's own
category label decides direction — case-insensitive substring match, so
`wc_gain`/`small_gain`/`mosaic_gain`/`amplification` all read as gain and
`wc_loss`/`(wc_)loss` as loss (`wc` = whole chromosome). A segment whose category
has none of those terms is ignored. This sidesteps ploidy entirely (a tetraploid
line's CN-2/3 losses are labelled `loss`, so they count as losses without any
baseline setting). `--cn-col` (default 13) is now recorded in the per-event table
only and never affects direction — a wrong/empty CN column can't zero out the
reference.

**inferCNV neutral state** depends on the HMM model: **i6** (states 1-6) is diploid
at **3**, **i3** (states 1-3) is diploid at **2**. The model is autodetected from
the query filename (`...predHMMi3...` vs `...HMMi6...`); override with
`--hmm-i {3,6}`. Getting this wrong flips gains/losses, so check the printed
`inferCNV model: HMMiN` line.

A reference event **matches** when same-direction inferCNV calls cover at least
`--min-overlap` (default **0.50**) of its length. Overlapping inferCNV intervals
are merged per direction first, so base pairs are never double-counted (this also
unions several `cell_group_name` groups cleanly; use `--group` to restrict to one,
or `plot_clones.py` to treat each as its own query).

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

## 3. Genome-wide overlay plot (one query)

This is the engine behind the per-clone overlays (step 4); run it directly for a
single clone (`--group <cell_group_name>`) or, without `--group`, for the union of
all clones. The all-samples pipeline does **not** emit a sample-level union overlay
(unioning subclones is coarse — the per-clone summary in step 4 is more useful).

```bash
python inferCNV/snp_array/plot_cnv_overlay.py \
    -r cellline_SKNBE2c.hg38.bed \
    -q results/BMO-SKNBE2c/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat \
    -o snp_vs_infercnv.overlay.png --group neuron.neuron_s1 \
    --title "BMO-SKNBE2c neuron_s1"
```

Chromosomes 1-22 laid end-to-end along x (add `--include-xy` for X/Y), with three
stacked bands:

- **reference** — SNP-array segments, gain = red, loss = blue.
- **query** — inferCNV HMM regions (union of all cell groups, or one via
  `--group`), same colours.
- **match** — each reference event decomposed into **matched** (green,
  same-direction inferCNV overlap), **contradicted** (purple, opposite call), and
  **missed** (grey, no inferCNV call). The title reports the overall % of
  reference bp matched / contradicted (identical to the `compare_cnv_overlap.py`
  "all" row).

`.png` or `.pdf` is chosen from the `-o` extension; `--dpi`, `--width`, `--height`
tune the raster. Large arm-level events dominate the picture; sub-Mb focal calls
are a few pixels wide (see resolution caveat below).

## 4. Per-clone plots (analysis_mode = "subclusters")

When inferCNV is run in `subclusters` mode, `pred_cnv_regions.dat` holds one
`cell_group_name` per clone. `plot_clones.py` treats **each clone as its own
query** and, in one call per sample, writes:

- one **3-band overlay per clone** (`<prefix>.<clone>.overlay.png`, like step 3) —
  skip with `--no-overlays`;
- a **clone summary** (`<prefix>.clones_summary.png`): every clone's *match* band
  stacked on shared genome axes (a reference band on top for context), rows sorted
  largest-clone-first, with a **clone-size sidebar**;
- a **table** (`<prefix>.clones.tsv`): per-clone cell count + % events / % bp
  matched / % bp contradicted.

```bash
python inferCNV/snp_array/plot_clones.py \
    -r cellline_SKNBE2c.hg38.bed \
    -q results/BMO-SKNBE2c/17_HMM_predHMMi3.leiden.hmm_mode-subclusters.pred_cnv_regions.dat \
    -g results/BMO-SKNBE2c/infercnv.17_HMM_predHMMi3.leiden.hmm_mode-subclusters.observation_groupings.txt \
    -o results/BMO-SKNBE2c/snp_vs_infercnv \
    --title BMO-SKNBE2c
```

Clone sizes come from `observation_groupings.txt`: its "Dendrogram Group" column
is the subcluster id (e.g. `neuron_s8`) and matches the suffix of a
`cell_group_name` (`neuron.neuron_s8`). Without `-g` the plots/table still work but
show no sizes (rows then sort by % bp matched). `--min-cells N` drops tiny clones.

## Running all samples — Snakemake pipeline (SLURM)

`Snakefile` runs the pipeline for every sample as a small standalone workflow that
reuses the repo's SLURM profile (`profiles/slurm`). Per sample it runs
`Liftover` -> `Compare` -> `Clones` (each its own SLURM job, in the
`envs/snp_array.yaml` conda env).

1. Edit `config_snp_array.yaml`: set `infercnv_results_dir`, `output_dir`, the
   `chain` path, and one `samples:` entry per sample pointing at its **hg19** array
   BED (samples that share a cell line — e.g. IMR DOX/noDOX — point at the same
   file; drop any sample you have no array for). A sample value may instead be a
   mapping `{bed:, type_col:, cn_col:}` to override the category / copy-number
   column for that one sample (rarely needed) — see the global `type_col` / `cn_col`
   and the comment above `samples:`.
2. Submit the controller from the `MultiomePipe/` root:

```bash
sbatch inferCNV/snp_array/run_snp_array_slurm.sh
# or a different config:
sbatch inferCNV/snp_array/run_snp_array_slurm.sh config/my_snp_array.yaml
```

Or run locally:

```bash
snakemake -s inferCNV/snp_array/Snakefile \
    --configfile inferCNV/snp_array/config_snp_array.yaml \
    --workflow-profile profiles/slurm --jobs 20
```

The controller needs a conda env with `snakemake` + `snakemake-executor-plugin-slurm`
(`CONDA_ENV` in the runner, default `snakemake`); each rule's own tools come from
`envs/snp_array.yaml`. Per sample the pipeline writes, under
`<output_dir>/<sample>/`: `snp_vs_infercnv.summary.tsv` / `.per_event.tsv`,
`.clones_summary.png` / `.clones.tsv`, and per-clone `.<clone>.overlay.png`. The
inferCNV pred / groupings files are found by glob (see `pred_glob` /
`groupings_glob`), preferring the Bayes-filtered `Pnorm` predictions.

To run one sample by hand instead, call the scripts directly as in sections 1-4
above.

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
- **Clones are inferCNV's own subclusters**, not independent truth. Small clones
  have noisy CNV profiles, so read the summary top-down (largest first) and weight
  by the size sidebar; a bulk SNP array can't confirm a rare subclone anyway.
- **Build.** Both files must be hg38 before step 2; don't skip the liftover.
