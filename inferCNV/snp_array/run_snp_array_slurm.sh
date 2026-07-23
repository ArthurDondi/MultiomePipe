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
# Submit from anywhere — the script cd's to its own inferCNV/snp_array/ dir so the
# relative -s Snakefile / --configfile config_snp_array.yaml resolve, and points the
# profile at the repo-root profiles/slurm:
#   sbatch inferCNV/snp_array/run_snp_array_slurm.sh
#   # or with a different config (relative to inferCNV/snp_array/, or absolute):
#   sbatch inferCNV/snp_array/run_snp_array_slurm.sh /abs/path/my_snp_array.yaml
#
# Needs a conda env with snakemake + the SLURM executor plugin (CONDA_ENV, default
# MultiomePipe). The plot/liftover tools (matplotlib, pyliftover) come from the
# per-rule env envs/snp_array.yaml, which snakemake builds because the profile sets
# use-conda: true — so matplotlib does NOT need to be in CONDA_ENV.

set -euo pipefail

# Resolve paths from the script's own location so it works from any submit dir.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$SCRIPT_DIR"          # so `-s Snakefile` and `--configfile config_snp_array.yaml` resolve

CONFIG="${1:-config_snp_array.yaml}"

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
echo "job id     : ${SLURM_JOB_ID:-<interactive>}"
echo "config     : $CONFIG"
echo "snakemake  : $(command -v snakemake)"
echo "======================"

# --workflow-profile profiles/slurm supplies the executor, per-rule queues and
# use-conda; --jobs caps how many rule-jobs run at once.
snakemake \
    -s Snakefile \
    --configfile "$CONFIG" \
    --workflow-profile "$REPO_ROOT/profiles/slurm" \
    --jobs 20 \
    --rerun-triggers mtime params software-env \
    -p
