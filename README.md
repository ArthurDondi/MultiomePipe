# MultiomePipe

MultiomePipe is a Snakemake workflow to analyze 10X Multiome data

## Methods

MultiomePipe is designed to process 10X Multiome ATAC + Gene Expression data after running Cellranger ARC. It also works with RNA-only or ATAC-only data (without integration).

### RNA-seq
MultiomePipe uses [Scanpy](https://scanpy.readthedocs.io/en/stable/index.html) to process, cluster, and annotate RNA-seq data. It was inspired by the [nf-core workflow scdownstream](https://github.com/nf-core/scdownstream),  the [single-cell best practices book](https://github.com/theislab/single-cell-best-practices), and [scanpy tutorials](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html).

Current steps:
- Corrects for ambient RNA and identify empty droplets using [CellBender](https://github.com/broadinstitute/CellBender)
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

### RNA-seq only
```
conda create -n MultiomePipe -c conda-forge -c bioconda snakemake python=3.11.4
conda activate MultiomePipe
git clone https://github.com/ArthurDondi/MultiomePipe.git
cd MultiomePipe/
```

### RNA-seq and ATAC-seq

If you also want to process ATAC-seq data, note that SCENIC+ could not be exported as a standalone conda environment for this pipeline. Therefore, you will need to install SCENIC+ directly into the main conda environment.

From MultiomePipe directory, with MultiomePipe conda env activated:
```
conda install -c bioconda pybedtools -y
git clone https://github.com/aertslab/scenicplus.git
cd scenicplus
```
In the scenicplus repo, remove the pybedtools lines in `requirements.in` and `requirements.txt`.

Change `requires-python = ">=3.8,<=3.11.11"` to `requires-python = ">=3.8,<=3.13.5"` in pyproject.toml (see [scenicplus#558](https://github.com/aertslab/scenicplus/pull/558))

Then:
```
pip install . --no-deps
pip install -r requirements.txt --no-deps
pip install git+https://github.com/macs3-project/MACS.git
pip install diptest
```

## Example with GEX only `config/config_Fetahu2023`
## Prepare

In `config/config_Fetahu2023`, change the `User` fields to paths to your desired input and output folders. 

Input folder needs to exist and to have a `marker_genes.json` file in it (to plot markers before manual annotation). For this example, you can copy it from `data/marker_genes.json`

## Run

```
snakemake -s workflow/Snakefile --configfile config/config_Fetahu2023.yaml --cores 1 --use-conda -p --resources mem_mb=16000 --conda-frontend conda -n
```

Note that for the moment, the pipeline can only run on a single core due to a bug with Cellbender temporary data overwriting (and to my limited current setup to test a fix).

Environments are all saved as `envs/***.yaml` files and will install automatically. You should be able to run the pipeline at least until requiring manual cluster annotations.



