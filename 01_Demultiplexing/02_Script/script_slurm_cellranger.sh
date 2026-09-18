#!/bin/bash
#SBATCH --job-name=cellranger_multi
#SBATCH --output=cellranger_multi_%j.out
#SBATCH --error=cellranger_multi_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --partition=debug
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lorena.quatreville@univ-amu.fr

set -euo pipefail

# ================================
# Container
# ================================
IMAGE=/data/madbio/image.sif
BASE_DIR=/data/madbio/Analysis
FASTQ_DIR=$BASE_DIR/Results/FASTQ
REF=$BASE_DIR/References/refdata-gex-GRCh38-2020-A

################ CellRanger

# Check if container exists
if [ ! -f "$IMAGE" ]; then
    echo "Image SIF introuvable : $IMAGE"
    exit 1
fi

# Map of runs to config CSVs
declare -A RUNS
RUNS[lib1]="$BASE_DIR/References/path_to_csv_file/Lib1/config.csv"
RUNS[lib2]="$BASE_DIR/References/path_to_csv_file/Lib2/config.csv"
RUNS[lib3]="$BASE_DIR/References/path_to_csv_file/Lib3/config.csv"
RUNS[lib4]="$BASE_DIR/References/path_to_csv_file/Lib4/config.csv"


LOCALMEM=120

for RUN in lib1 lib2 lib3 lib4; do
    CONFIG="${RUNS[$RUN]}"
    OUTDIR="$BASE_DIR/Results/cellranger_multi"
    mkdir -p "$OUTDIR"

    if [ ! -f "$CONFIG" ]; then
        echo "Fichier config introuvable : $CONFIG"
        exit 1
    fi

    echo "=== Lancement de $RUN ==="
    echo "Config : $CONFIG"
    echo "Outdir : $OUTDIR"
    echo "Cores  : $SLURM_CPUS_PER_TASK"
    echo "Mem    : ${LOCALMEM}G"

    cd "$OUTDIR"

    singularity exec \
        --bind $BASE_DIR:$BASE_DIR \
        $IMAGE \
        cellranger multi \
        --id="$RUN" \
        --csv="$CONFIG" \
        --localcores="$SLURM_CPUS_PER_TASK" \
        --localmem="$LOCALMEM"

    echo "$RUN terminé avec succès."

    cd - >/dev/null
done

echo "Tous les runs ont été traités avec succès."
