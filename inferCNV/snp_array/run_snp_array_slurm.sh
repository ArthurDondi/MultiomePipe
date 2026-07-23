#!/bin/bash
#SBATCH --output /nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/logs/snp_array_%j.log
#SBATCH --error  /nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/logs/snp_array_%j.err
#SBATCH --job-name=snp_array
#SBATCH --partition=shortq     # 12h is ample: this controller only submits a handful of quick jobs
#SBATCH --qos=shortq           # qos must match the partition
#SBATCH --time=0-06:00:00      # --time=days-hours:minutes:seconds
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      # this job is only the Snakemake controller
#SBATCH --mem=4000             # controller is light; rules get their own jobs
#SBATCH --mail-type=end
#SBATCH --mail-user=arthur.dondi@cemm.at
#
# Run the SNP-array vs inferCNV comparison pipeline on SLURM.
#
# This batch job is the Snakemake *controller*: it stays alive for the run and
# submits each rule (Liftover / Compare / Overlay / Clones, per sample) as its own
# SLURM job through profiles/slurm/config.yaml. It therefore needs few resources
# itself. Each rule uses the inferCNV/snp_array/envs/snp_array.yaml conda env.
#
# Submit from inside the repo (typically the MultiomePipe/ root, like
# run_Multiomepipe_slurm.sh); the script locates the repo root itself:
#   sbatch inferCNV/snp_array/run_snp_array_slurm.sh
#   # or with a different config (absolute, repo-root-relative, or a bare filename):
#   sbatch inferCNV/snp_array/run_snp_array_slurm.sh /abs/path/my_snp_array.yaml
#
# Needs a conda env with snakemake + the SLURM executor plugin (CONDA_ENV, default
# MultiomePipe). The plot/liftover tools (matplotlib, pyliftover) come from the
# per-rule env envs/snp_array.yaml, which snakemake builds because the profile sets
# use-conda: true — so matplotlib does NOT need to be in CONDA_ENV.

set -euo pipefail

# Find the repo root (the dir containing profiles/slurm) by walking up from the
# submit dir. NB: do NOT use $BASH_SOURCE here — under sbatch, SLURM copies this
# script to a spool dir, so the script's own path is not inside the repo.
start_dir="${SLURM_SUBMIT_DIR:-$PWD}"
REPO_ROOT="$start_dir"
while [ "$REPO_ROOT" != "/" ] && [ ! -d "$REPO_ROOT/profiles/slurm" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
if [ ! -d "$REPO_ROOT/profiles/slurm" ]; then
    echo "ERROR: could not find the MultiomePipe root (a dir with profiles/slurm)"
    echo "       by walking up from '$start_dir'. Submit from inside the repo, e.g.:"
    echo "         cd /path/to/MultiomePipe && sbatch inferCNV/snp_array/run_snp_array_slurm.sh"
    exit 1
fi
cd "$REPO_ROOT"

# Config: default, or $1 given as an absolute path / repo-root-relative / bare name.
CONFIG="${1:-inferCNV/snp_array/config_snp_array.yaml}"
if [ ! -f "$CONFIG" ] && [ -f "inferCNV/snp_array/$CONFIG" ]; then
    CONFIG="inferCNV/snp_array/$CONFIG"
fi
if [ ! -f "$CONFIG" ]; then
    echo "ERROR: config file not found: '$CONFIG' (from repo root $REPO_ROOT)"; exit 1
fi

# --- activate the conda env that has snakemake + the slurm executor plugin ----
CONDA_BASE="${CONDA_BASE:-/nobackup/lab_taschner-mandl/arthurdondi/miniconda3}"
CONDA_ENV="${CONDA_ENV:-MultiomePipe}"
# conda's shell hook is not always `set -u` clean, so relax nounset around it.
set +u
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"
set -u

echo "======================"
echo "submit dir : ${SLURM_SUBMIT_DIR:-$(pwd)}"
echo "repo root  : $REPO_ROOT"
echo "job id     : ${SLURM_JOB_ID:-<interactive>}"
echo "config     : $CONFIG"
echo "snakemake  : $(command -v snakemake)"
echo "======================"

# Everything below is repo-root-relative (we cd'd there). --workflow-profile
# profiles/slurm supplies the executor, per-rule queues and use-conda; --jobs caps
# how many rule-jobs run at once.
snakemake \
    -s inferCNV/snp_array/Snakefile \
    --configfile "$CONFIG" \
    --workflow-profile profiles/slurm \
    --jobs 20 \
    --rerun-triggers mtime params software-env \
    -p
