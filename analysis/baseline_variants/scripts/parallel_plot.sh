#!/usr/bin/env bash

pyscript="$1"
indir="$2"
threads="$3"

fd -j "$threads" -HI -t d -d 1 \
    -x python "$pyscript" '{}' '{}' \; \
    . "$indir"