MultiomePipe is a Snakemake workflow using scverse (scanpy) to pre-process 10X Multiome data after running Cellranger ARC. Alternatively, it also works with GEX-only data.

This workflow was inspired by the [nf-core workflow scdownstream](https://github.com/nf-core/scdownstream),  the [single-cell best practices book](https://github.com/theislab/single-cell-best-practices), and scanpy [tutorials](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html).


# Install
```
conda create -n MultiomePipe -c conda-forge -c bioconda snakemake python=3.11.4
conda activate MultiomePipe
git clone https://github.com/ArthurDondi/MultiomePipe.git
cd MultiomePipe/
```

# Example with GEX only `config/config_Fetahu2023`
## Prepare

In `config/config_Fetahu2023`, change the `User` fields to paths to your desired input and output folders. 

Input folder needs to exist and to have a `marker_genes.json` file in it (to plot markers before manual annotation). For this example, you can copy it from `data/marker_genes.json`

## Run

```
snakemake -s workflow/Snakefile --configfile config/config_Fetahu2023.yaml --cores 1 --use-conda -p --resources mem_mb=16000 --conda-frontend conda -n
```

Note that for the moment, the pipeline can only run on a single core due to a bug with Cellbender temporary data overwriting (and to my limited current setup to test a fix).

Environments are all saved as `envs/***.yaml` files and will install automatically. You should be able to run the pipeline at least until requiring manual cluster annotations.



