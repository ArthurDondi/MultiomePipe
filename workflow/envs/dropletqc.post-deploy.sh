#!/usr/bin/env bash
# Snakemake post-deploy script for envs/dropletqc.yaml.
#
# Runs once, single-threaded, immediately after Snakemake creates the conda
# environment (with the env activated and $CONDA_PREFIX set). It installs the
# DropletQC package from GitHub into the env's R library so that:
#   * the install never races against parallel DropletQC jobs (no 00LOCK errors),
#   * every dependency is already present as a conda binary (no source builds),
#   * the package lives on the env's default .libPaths(), so the multisession
#     future workers spawned by nuclear_fraction_tags() can attach it.
#
# All runtime dependencies are provided by dropletqc.yaml, so we install with
# dependencies = FALSE and upgrade = "never" for a deterministic, offline-safe
# install.
set -euo pipefail

Rscript -e '
  if (!requireNamespace("DropletQC", quietly = TRUE)) {
    remotes::install_github(
      "powellgenomicslab/DropletQC",
      upgrade      = "never",
      dependencies = FALSE
    )
  }
  # Fail the env build loudly if the package still cannot be loaded.
  library(DropletQC)
  cat("DropletQC installed:", as.character(packageVersion("DropletQC")), "\n")
'
