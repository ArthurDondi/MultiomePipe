#!/bin/bash
# Run MultiomePipe on a SLURM cluster.
#
# Each rule is submitted as an individual Slurm job. Default resources
# (partition, memory, runtime, CPUs) are defined in profile/slurm/config.yaml.
# Edit that file to match your cluster's configuration before running.
#
# Requirements:
#   pip install snakemake-executor-plugin-slurm
#
# Usage (from the MultiomePipe/ root directory):
#   bash run_Multiomepipe_slurm.sh
#
# Slurm logs are written to logs/slurm/ inside the output directory.
#
# --configfile                 : change to point to your own config file
# --profile profile/slurm      : load the Slurm executor profile
# --resources nvidia_gpu=1     : declare 1 GPU available; the
#                                CellbenderRemoveBackgroundRNA rule requests
#                                nvidia_gpu=1, so at most one CellBender job
#                                runs at a time.
# --rerun-triggers mtime params: resubmit when inputs or params change


snakemake \
    -s workflow/Snakefile \
    --configfile config/config_BMO_Viki.yaml \
    --profile profile/slurm \
    --resources nvidia_gpu=1 \
    --rerun-triggers mtime params \
    -p
