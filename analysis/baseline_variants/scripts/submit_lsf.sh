#!/usr/bin/env bash
set -eux

JOB_NAME="snakemake_master_process."$(date --iso-8601='minutes')
LOG_DIR=logs/

if [[ ! -d "$LOG_DIR" ]]; then
    echo "Error: Log directory $LOG_DIR does not exist"
    exit 1
fi

MEMORY=4000
THREADS=4
PROFILE="lsf-dev"
BINDS="/tmp,$HOME"
case $HOSTNAME in
    *noah*)
        BINDS+=",/scratch,/hps/nobackup/research/zi,/nfs/research1/zi"
        ;;
    *codon*)
        BINDS+=",/hps/scratch,/hps/nobackup/iqbal,/nfs/research/zi,$FASTSW_DIR --scratch /hps/scratch"
        ;;
    *)
        echo "ERROR: HOSTNAME $HOSTNAME not recognised"
        exit 1
        ;;
esac

ARGS="--contain -B $BINDS"

bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY] span[hosts=1]" \
    -M "$MEMORY" \
    -n "$THREADS" \
    -o "$LOG_DIR"/"$JOB_NAME".o \
    -e "$LOG_DIR"/"$JOB_NAME".e \
    -J "$JOB_NAME" \
    snakemake --profile "$PROFILE" \
    --local-cores "$THREADS" \
    "$@" --singularity-args "$ARGS"

exit 0
