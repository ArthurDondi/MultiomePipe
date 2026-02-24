#!/bin/bash
# Run MultiomePipe pipeline in the background.
#
# --cores 8                 : up to 8 CPU cores in total
# --resources nvidia_gpu=1  : declare 1 GPU available; the CellbenderRemoveBackgroundRNA
#                             rule requests nvidia_gpu=1, so at most one CellBender job runs
#                             at a time while all other rules can still run in parallel.
#
# nohup ... &> snakemake.log & : run detached from the terminal; both stdout and stderr
#                               are written to snakemake.log; the job survives terminal close.
# $! is the PID of the background process, printed so you can monitor or kill it.

nohup snakemake \
    -s workflow/Snakefile \
    --configfile config/config_BMO_Viki.yaml \
    --cores 8 \
    --use-conda \
    -p \
    --resources mem_mb=64000 nvidia_gpu=1 \
    --conda-frontend conda \
    &> snakemake.log &

echo "Snakemake started in background (PID $!). Follow progress with: tail -f snakemake.log"
