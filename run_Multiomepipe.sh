#!/bin/bash
# Run MultiomePipe pipeline locally.
# --cores 8            : up to 8 CPU cores in total
# --jobs 1             : run only one job at a time
# --resources nvidia_gpu=1 : declare 1 GPU available; the CellbenderRemoveBackgroundRNA
#                       rule requests nvidia_gpu=1, so only one cellbender job will run
#                       at a time and it will be the only job running (combined with --jobs 1)

snakemake \
    -s workflow/Snakefile \
    --configfile config/config_BMO_Viki.yaml \
    --cores 8 \
    --jobs 1 \
    --use-conda \
    -p \
    --resources mem_mb=64000 nvidia_gpu=1 \
    --conda-frontend conda

exit $?
