# MultiomePipe

A Snakemake workflow to analyze 10X Multiome data

## Methods

MultiomePipe is designed to process 10X Multiome ATAC + Gene Expression data after running Cellranger ARC. It also works with RNA-only or ATAC-only data (without integration).

The workflow has only been tested on human data and will require changes for other species.

### RNA-seq
MultiomePipe uses [Scanpy](https://scanpy.readthedocs.io/en/stable/index.html) to process, cluster, and annotate RNA-seq data. It was inspired by the [nf-core workflow scdownstream](https://github.com/nf-core/scdownstream),  the [single-cell best practices book](https://github.com/theislab/single-cell-best-practices), and [scanpy tutorials](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html).

Current steps:
- Corrects ambient RNA either with [CellBender](https://github.com/broadinstitute/CellBender) or with DropletQC damaged-cell filtering followed by [SoupX](https://github.com/constantAmateur/SoupX)
- Filter cells with fewer than (default 200) genes expressed
- Filter out cells whose mitochondrial content deviates from the median by more than 5 MADs
- Identify and remove doublets using [Scrublet](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.scrublet.html)
- Normalize counts with median count depth normalization
- Correct for batches using [Harmony](https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pp.harmony_integrate.html)
- Embedding using UMAP and clustering using Leiden
- Plots user-provided marker genes for manual cell annotation

### ATAC-seq
The ATAC-seq processing in this workflow is an implementation of [SCENIC+](https://www.nature.com/articles/s41592-023-01938-4) (only pyCisTopic for now).

The pyCisTopic implementation largely follows this [pyCisTopic Tutorial](https://pycistopic.readthedocs.io/en/latest/notebooks/human_cerebellum.html), with MACS3 instead of MACS2.

## Install
Create a conda environment
```
conda create -n MultiomePipe python=3.11 -y
conda activate MultiomePipe
```
Clone the MultiomePipe git
```
git clone https://github.com/ArthurDondi/MultiomePipe.git
cd MultiomePipe/
```
Note that SCENIC+ could not be exported as a standalone conda environment for this pipeline. Therefore, you will need to install SCENIC+ directly into the main conda environment.

Clone the SCENIC+ git (doesn't need to be in the MultiomPipe git)
```
git clone https://github.com/aertslab/scenicplus.git
cd scenicplus
```
In the scenicplus dir, remove the pybedtools lines in `requirements.in` and `requirements.txt`.

Change `requires-python = ">=3.8,<=3.11.11"` to `requires-python = ">=3.8,<=3.13.5"` in pyproject.toml (see [scenicplus#558](https://github.com/aertslab/scenicplus/pull/558))

Now, with the MultiomePipe conda env activated and from the scenicplus directory:
```
conda install -c bioconda pybedtools -y
pip install . --no-deps
pip install -r requirements.txt --no-deps
pip install git+https://github.com/macs3-project/MACS.git
pip install diptest
```

## Examples
### ATAC-seq + RNA-seq: `config/config_pbmc_unsorted_3k`
#### Prepare

In `config/config_pbmc_unsorted_3k`, change the `User` fields to paths to your desired input and output folders. 

Input folder needs to have a `marker_genes.json` (exact name) file in it, to plot markers before manual annotation. For this example, you can copy it from `MultiomePipe/data/marker_genes_pbmc_unsorted_3k.json`.

#### Run
The pipeline should be called from `MultiomePipe/`. First, do a dry run (`-n`):
```
snakemake -s workflow/Snakefile --configfile config/config_pbmc_unsorted_3k.yaml --cores 1 --use-conda -p --resources mem_mb=16000 --conda-frontend conda -n
```

If it works, you can run:
```
snakemake -s workflow/Snakefile --configfile config/config_pbmc_unsorted_3k.yaml --cores 1 --use-conda -p --resources mem_mb=16000 --conda-frontend conda 
```

Note that for the moment, the pipeline can only run on a single core due to a bug with Cellbender temporary data overwriting (and to my limited current setup to test a fix).

Environments are all saved as `envs/***.yaml` files and will install automatically. You should be able to run the pipeline until requiring manual cluster annotations.

#### Manual annotation
The rule `CheckManualAnnotation` will fail until you fill the `youroutputfolder/QC/RNA/Merged/Annotation/merged_manual_annotation.json` file in the following format:
```
{
  "0": "T cell",
  "1": "Monocytes",
  "2": "T cell",
  "3": "Monocytes",
  "4": "B cell",
  "5": "T cell",
  "6": "NK cell",
  "7": "Myeloid cell",
  "8": "Dendritic cell",
  "9": "pDCs",
  "10": "Myeloid cell"
}
```
With `"leiden_cluster_number": "cell_type"`.

For this, you can use  the files `umap_leiden_res_merge.png` and `dotplot__manual_markers_leiden_res_{0.20/0.50/1.00}_merge.png` in `youroutputfolder/QC/RNA/Merged/BatchCorrection/`. You can select your favorite Leiden resolution (0.2, 0.5, or 1.00), but don't forget to match `leiden_res` under `PlottingAnnotations` in your `config.yaml`.

Once `youroutputfolder/QC/RNA/Merged/Annotation/merged_manual_annotation.json` is filled, you can re-run the pipeline with the same command.

### RNA-seq only: `config/config_Fetahu2023`
#### Prepare

In `config/config_Fetahu2023`, change the `User` fields to paths to your desired input and output folders. 

Input folder needs to have a `marker_genes.json` (exact name) file in it, to plot markers before manual annotation. For this example, you can copy it from `MultiomePipe/data/marker_genes_Fetahu2023.json`

#### Run
```
snakemake -s workflow/Snakefile --configfile config/config_Fetahu2023.yaml --cores 1 --use-conda -p --resources mem_mb=16000 --conda-frontend conda 
```

#### Manual annotation
See above.

## Running on a SLURM cluster

MultiomePipe can be submitted to a SLURM cluster using the profile in
`profiles/slurm/config.yaml` (Snakemake 8 executor plugin).

### One-time setup
Install the SLURM executor plugin into the MultiomePipe conda env:
```
pip install snakemake-executor-plugin-slurm
```

### Edit the profile for your cluster
Open `profiles/slurm/config.yaml` and adjust the partitions / QoS to match your
site. The defaults target the queues:
- `shortq` (12h) for the default of most rules,
- `mediumq` (2d) for `CellRangerCount` and `RunLDAModels`,
- `gpu` for `CellbenderRemoveBackgroundRNA` (requests `--gres=gpu:1`).

On this cluster the qos must match the partition, which is why each partition is
paired with its matching `qos` resource. QoS is set via the dedicated `qos`
resource (and GPUs via `gres`) rather than `slurm_extra`, because the SLURM
executor plugin manages those itself and rejects `--qos` / `--gres` passed
through `slurm_extra`. `runtime` is expressed in **minutes**. Per-rule `mem_mb` /
`runtime` are declared on the rules in `workflow/rules/*.smk`; the profile only
routes heavy rules to the right queue. Tune `mem_mb` / `runtime` per rule if your
data is larger.

### Submit
From `MultiomePipe/` (point input/output dirs in your config at a shared
filesystem such as `/nobackup/<group>/<user>/...`, not your home), first do a
dry run on the login node:
```
snakemake -s workflow/Snakefile \
    --configfile config/<your_config>.yaml \
    --workflow-profile profiles/slurm \
    -n          # dry run; drop -n to submit interactively
```

To run the whole workflow unattended, submit the controller batch job
`run_Multiomepipe_slurm.sh`, which runs Snakemake on `longq` (30d) and submits
every rule as its own SLURM job via the profile above:
```
sbatch run_Multiomepipe_slurm.sh config/<your_config>.yaml
```
Edit the `#SBATCH` log paths / `--mail-user` at the top of that script for your
account before submitting.

### Caveats
- **CellBender parallelism:** the `CellbenderRemoveBackgroundRNA` rule writes a
  `ckpt.tar.gz` in the shared working directory, so running multiple CellBender
  jobs concurrently can clobber each other. Until this is refactored to a
  per-job temp dir, keep CellBender effectively serial (e.g. submit one sample
  at a time, or set `--cores`/`jobs` low for that step).
- **Tools on compute nodes:** `cellranger` (path in `config['User']['cellranger']`)
  and Mallet must be available / `module load`-able on the compute nodes, not
  only the login node.