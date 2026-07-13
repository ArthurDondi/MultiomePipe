# inferCNV from a MultiomePipe batch-corrected object

Downstream (manual) copy-number inference on the RNA cohort, run in **R / RStudio**
after the `BatchCorrection` step. Not wired into the Snakemake workflow — it is a
two-step bridge: export the `.h5ad` to flat files (Python), then run inferCNV (R).

```
BatchCorrection  ->  merged.batch_corrected.h5ad
        |
        |  1) run_export_slurm.sh -> export_for_infercnv.py   (scverse env, SLURM)
        v
   inferCNV/input/ : counts.mtx.gz, genes.tsv, barcodes.tsv,
                     metadata.tsv, annotations_celltype.tsv, gene_ordering.tsv
        |
        |  2) run_infercnv.R           (RStudio)
        v
   inferCNV/results/<sample>/ : infercnv.png, HMM_CNV_predictions.*, ...
```

## 1. Export (Python, on the cluster)

Reads raw counts from `layers['counts']`, cell groups from `C_scANVI` (the best
label available at the BatchCorrection stage — there is no manual `cell_type`
yet), and builds the gene→chromosome ordering from the same GTF used for the
CellRanger reference. All three paths default to `config/config_BMO_combined.yaml`
values, so on the cluster you usually need no arguments.

**Submit as a SLURM job (recommended).** `run_export_slurm.sh` activates the
`scverse` env and runs the exporter on the `shortq` queue (2 cpus / 64 GB / 4 h):

```bash
# from the MultiomePipe/ root
sbatch inferCNV/run_export_slurm.sh
```

Any flags are forwarded straight to `export_for_infercnv.py`, so you can override
the defaults — e.g. a different output dir, or export the *annotated* object:

```bash
sbatch inferCNV/run_export_slurm.sh --outdir /nobackup/.../inferCNV/input_v2
sbatch inferCNV/run_export_slurm.sh --input /path/merged.annotated.h5ad --annotation-key cell_type
```

Conda installed elsewhere? `CONDA_BASE=/path/to/miniconda3 CONDA_ENV=scverse sbatch ...`.
Logs land in `projects/BMO/logs/inferCNV_export_<jobid>.log`; the exporter prints
the chosen label column and the cells-per-group table — check that the reference
cell types in `run_infercnv.R` match.

**Or run it directly** (e.g. on an interactive node), same flags:

```bash
conda activate scverse
python inferCNV/export_for_infercnv.py            # uses the config defaults
```

## 2. Run inferCNV (R / RStudio)

One-time install:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("infercnv")     # + Matrix
# HMM step also needs JAGS: https://mcmc-jags.sourceforge.io/  (conda: conda install -c conda-forge jags)
```

Edit the `USER SETTINGS` block at the top of `run_infercnv.R` (`INPUT_DIR`,
`OUTPUT_DIR`, and `NORMAL_CELL_TYPES` — the diploid reference), then run it. With
**28 samples** the default `MODE = "per_sample"` loops over each sample and uses
that sample's own normal (T / NK / B / myeloid / erythroid / pDC / HSPC) cells as
the reference; malignant *Neuroblastoma cell*s become the observations. Pure
cell-line samples with no normal cells fall back to a reference-less run
(baseline = mean of all cells). Set `MODE = "joint"` to run the whole cohort at
once.

Key knobs: `ANALYSIS_MODE` (`"samples"` fast / `"subclusters"` for tumour
subclones), `USE_HMM`, `NUM_THREADS`, and `SUBSAMPLE_MAX` (cap cells per group
for speed).

## Notes

- inferCNV needs **raw** counts (it models expression itself) — the exporter uses
  `layers['counts']`, not the log-normalised `X`.
- The matrix is exported sparse (gzipped MatrixMarket, genes × cells) so the full
  cohort stays small; `run_infercnv.R` reads it with `Matrix::readMM()`.
- The GFP transgene and any scaffold/`chrM` genes have no standard chromosomal
  position, so they are left out of `gene_ordering.tsv` and inferCNV ignores them.
- `CreateInfercnvObject` excludes chrX/chrY/chrM by default; edit `chr_exclude`
  in `run_infercnv.R` if you want the sex chromosomes.
```
