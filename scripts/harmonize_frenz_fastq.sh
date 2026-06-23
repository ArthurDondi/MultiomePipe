#!/bin/bash
#SBATCH --output /nobackup/lab_taschner-mandl/arthurdondi/harmonize_frenz_%A_%a.log
#SBATCH --error  /nobackup/lab_taschner-mandl/arthurdondi/harmonize_frenz_%A_%a.err
#SBATCH --job-name=harmonize_frenz
#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=8000
#SBATCH --array=1-10

set -euo pipefail

SRA_DIR=/home/arthur_d/bioinf_isilon/Research/TASCHNERMANDL/Internal/data/frenz-wiessner2024.scRNA/fastq
OUT_BASE=/path/to/input_fastq/Frenz    # <-- change this

# SRR  sample  lane
MAP=(
  "SRR26993577 Frenz-Wiessner_BMO5 1"
  "SRR26993578 Frenz-Wiessner_BMO5 2"
  "SRR26993579 Frenz-Wiessner_BMO4 1"
  "SRR26993580 Frenz-Wiessner_BMO4 2"
  "SRR26993581 Frenz-Wiessner_BMO3 1"
  "SRR26993582 Frenz-Wiessner_BMO3 2"
  "SRR26993583 Frenz-Wiessner_BMO2 1"
  "SRR26993584 Frenz-Wiessner_BMO2 2"
  "SRR26993585 Frenz-Wiessner_BMO1 1"
  "SRR26993586 Frenz-Wiessner_BMO1 2"
)

read -r SRR SAMPLE LANE <<< "${MAP[$((SLURM_ARRAY_TASK_ID-1))]}"
LANE_TAG=$(printf 'L%03d' "$LANE")
OUTDIR="$OUT_BASE/$SAMPLE"
mkdir -p "$OUTDIR"

echo "[$(date)] SRR=$SRR  sample=$SAMPLE  lane=$LANE_TAG"

# _1 = R1 (16bp CB + 12bp UMI) ; _2 = R2 (cDNA)
pigz -p "${SLURM_CPUS_PER_TASK}" -c "$SRA_DIR/${SRR}_1.fastq" > "$OUTDIR/${SAMPLE}_S1_${LANE_TAG}_R1_001.fastq.gz"
pigz -p "${SLURM_CPUS_PER_TASK}" -c "$SRA_DIR/${SRR}_2.fastq" > "$OUTDIR/${SAMPLE}_S1_${LANE_TAG}_R2_001.fastq.gz"

echo "[$(date)] done -> $OUTDIR"
