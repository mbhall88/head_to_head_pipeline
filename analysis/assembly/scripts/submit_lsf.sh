#!/usr/bin/env bash
JOB_NAME="snakemake_master_process."$(date --iso-8601='minutes')
LOG_DIR=logs/

if [[ ! -d "$LOG_DIR" ]]
then
    echo "Error: Log directory $LOG_DIR does not exist"
    exit 1
fi

MEMORY=4000
PROFILE="lsf"

#lustre_mount="-B /hps/nobackup/research/zi/"
module load singularity/3.5.0

bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY] span[hosts=1]" \
    -M "$MEMORY" \
    -o "$LOG_DIR"/"$JOB_NAME".o \
    -e "$LOG_DIR"/"$JOB_NAME".e \
    -J "$JOB_NAME" \
        snakemake --profile "$PROFILE" "$@"
#            --singularity-args "$lustre_mount" "$@"

exit 0