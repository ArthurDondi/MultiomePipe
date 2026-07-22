#!/bin/bash
#SBATCH --output /nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/logs/inferCNV_export_%j.log
#SBATCH --error  /nobackup/lab_taschner-mandl/arthurdondi/projects/BMO/logs/inferCNV_export_%j.err
#SBATCH --job-name=inferCNV_export
#SBATCH --partition=shortq      # 12h is plenty: this only reads the h5ad and writes flat files
#SBATCH --qos=shortq            # qos must match the partition on this cluster
#SBATCH --time=0-04:00:00       # --time=days-hours:minutes:seconds
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2       # exporter is single-threaded; a little headroom
#SBATCH --mem=64000             # holds the merged object + the transposed counts copy
#SBATCH --mail-type=end
#SBATCH --mail-user=arthur.dondi@cemm.at
#
# Export the batch-corrected .h5ad into inferCNV inputs on the SLURM cluster by
# running inferCNV/export_for_infercnv.py in the `scverse` conda env (which has
# anndata / scipy / pandas). Use this instead of running the exporter by hand.
#
# Submit from the MultiomePipe/ root:
#   sbatch inferCNV/run_export_slurm.sh
#   # or override any exporter flag (defaults come from config_BMO_combined.yaml):
#   sbatch inferCNV/run_export_slurm.sh --outdir /nobackup/.../inferCNV/input_v2
#   sbatch inferCNV/run_export_slurm.sh --input /path/merged.annotated.h5ad --annotation-key cell_type
#
# Conda install lives elsewhere? Override at submit time:
#   CONDA_BASE=/path/to/miniconda3 CONDA_ENV=scverse sbatch inferCNV/run_export_slurm.sh

set -euo pipefail

# --- activate the conda env that has anndata / scipy / scanpy ----------------
CONDA_BASE="${CONDA_BASE:-/nobackup/lab_taschner-mandl/arthurdondi/miniconda3}"
CONDA_ENV="${CONDA_ENV:-scverse}"
# conda's shell hook is not always `set -u` clean, so relax nounset around it.
set +u
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"
set -u

# Resolve export_for_infercnv.py next to this script, so cwd does not matter.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "======================"
echo "job id    : ${SLURM_JOB_ID:-<interactive>}"
echo "node      : $(hostname)"
echo "conda env : $CONDA_ENV ($CONDA_BASE)"
echo "python    : $(command -v python)"
echo "exporter  : $SCRIPT_DIR/export_for_infercnv.py"
echo "extra args: $*"
echo "======================"

python -u -W ignore export_for_infercnv.py "$@"
