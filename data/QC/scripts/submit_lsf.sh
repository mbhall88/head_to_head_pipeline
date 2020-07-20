#!/usr/bin/env bash
set -eux

JOB_NAME="snakemake_master_process."$(date --iso-8601='minutes')
LOG_DIR=logs/

if [[ ! -d "$LOG_DIR" ]]
then
    echo "Error: Log directory $LOG_DIR does not exist"
    exit 1
fi

MEMORY=4000
THREADS=4
PROFILE="lsf"
SINGULARITY_BINDS="/hps/nobackup/research/zi,/nfs/research1/zi"
SINGULARITY_WORKDIR="/scratch"
SINGULARITY_ARGS="--contain --workdir $SINGULARITY_WORKDIR --bind $SINGULARITY_BINDS --pwd $(pwd)"

echo "Passing the following args to singularity: $SINGULARITY_ARGS"

module load singularity/3.5.0

bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY] span[hosts=1]" \
    -M "$MEMORY" \
    -n "$THREADS" \
    -o "$LOG_DIR"/"$JOB_NAME".o \
    -e "$LOG_DIR"/"$JOB_NAME".e \
    -J "$JOB_NAME" \
snakemake --verbose --profile "$PROFILE" \
    --local-cores "$THREADS" \
    "$@" \
    --singularity-args "$SINGULARITY_ARGS"

exit 0