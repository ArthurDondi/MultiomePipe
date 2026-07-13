# inferCNV from a MultiomePipe batch-corrected object

Downstream (manual) copy-number inference on the RNA cohort, run in **R / RStudio**
after the `BatchCorrection` step. Not wired into the Snakemake workflow — it is a
two-step bridge: export the `.h5ad` to flat files (Python), then run inferCNV (R).

```
BatchCorrection  ->  merged.batch_corrected.h5ad
        |
        |  1) export_for_infercnv.py   (scverse conda env, on the cluster)
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
CellRanger reference. Defaults are filled from `config/config_BMO_combined.yaml`:

```bash
conda activate scverse
python inferCNV/export_for_infercnv.py \
    --input  /nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/QC/RNA/Merged/BatchCorrection/merged.batch_corrected.h5ad \
    --outdir /nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/inferCNV/input \
    --gtf    /nobackup/lab_taschner-mandl/arthurdondi/resources/references/hg38/BMO/gencode.v50.basic.annotation.plusGFP.gtf
```

Useful flags: `--annotation-key cell_type` (if you export the *annotated* object
instead), `--counts-layer counts`. The script prints the chosen label column and
the cells-per-group table — check that the reference cell types below match.

### From JupyterLab

Don't paste the script body into a cell and rely on the CLI — argparse would try
to parse the kernel's own `-f <kernel>.json` and fail. Either run it as a script
(`sys.argv` is set correctly, defaults come from the config):

```python
%run /path/to/MultiomePipe/inferCNV/export_for_infercnv.py            # uses config defaults
# %run .../export_for_infercnv.py --input ... --outdir ... --gtf ...  # or explicit paths
```

or import the function and call it with keyword arguments:

```python
import sys; sys.path.insert(0, "/path/to/MultiomePipe/inferCNV")
from export_for_infercnv import run_export
run_export(
    input_file="/nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/QC/RNA/Merged/BatchCorrection/merged.batch_corrected.h5ad",
    outdir="/nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/inferCNV/input",
    gtf="/nobackup/lab_taschner-mandl/arthurdondi/resources/references/hg38/BMO/gencode.v50.basic.annotation.plusGFP.gtf",
)
```

(The CLI now also ignores Jupyter's injected `-f` if you do paste-and-run.)

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
