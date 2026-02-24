#!/bin/bash
# Run MultiomePipe pipeline locally.
# --cores 8                 : up to 8 CPU cores in total
# --resources nvidia_gpu=1  : declare 1 GPU available; the CellbenderRemoveBackgroundRNA
#                             rule requests nvidia_gpu=1, so at most one CellBender job runs
#                             at a time while all other rules can still run in parallel.

snakemake \
    -s workflow/Snakefile \
    --configfile config/config_BMO_Viki.yaml \
    --cores 8 \
    --use-conda \
    -p \
    --resources mem_mb=64000 nvidia_gpu=1 \
    --conda-frontend conda

exit $?
