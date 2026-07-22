#!/bin/bash
#SBATCH --output /nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/logs/MultiomePipe_%j.log
#SBATCH --error  /nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/logs/MultiomePipe_%j.err
#SBATCH --job-name=MultiomePipe 
#SBATCH --partition=longq      # 30d limit: long enough for the whole workflow
#SBATCH --qos=longq            # qos must match the partition
#SBATCH --time=1-00:00:00     # --time=days-hours:minutes:seconds 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      # this job is only the Snakemake controller
#SBATCH --mem=8000             # controller is light; rules get their own jobs
#SBATCH --mail-type=end
#SBATCH --mail-user=arthur.dondi@cemm.at
#
# Run the whole MultiomePipe Snakemake workflow on the SLURM cluster.
#
# This batch job is the Snakemake *controller*: it stays alive for the whole
# run and submits each rule as its own SLURM job through the profile in
# profiles/slurm/config.yaml (per-rule memory / walltime / partition, and the
# GPU partition for CellBender). It therefore needs few resources itself but a
# long walltime, hence longq.
#
# Submit with:
#   sbatch run_Multiomepipe_slurm.sh [config/your_config.yaml]
#
# Must be submitted from the MultiomePipe/ root, with the MultiomePipe conda env
# activated (and `pip install snakemake-executor-plugin-slurm` done once).

set -euo pipefail

CONFIG="${1:-config/config_BMO_combined.yaml}"

echo "======================"
echo "submit dir : $SLURM_SUBMIT_DIR"
echo "job name   : $SLURM_JOB_NAME"
echo "partition  : $SLURM_JOB_PARTITION"
echo "job id     : $SLURM_JOB_ID"
echo "config     : $CONFIG"
echo "======================"

snakemake \
    -s workflow/Snakefile \
    --configfile "$CONFIG" \
    --workflow-profile profiles/slurm \
    --jobs 50 \
    --rerun-triggers mtime params software-env code \
    -p --until PlottingAnnotationsManual -R BatchCorrection
