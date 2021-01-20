#!/usr/bin/env bash
set -eu

infile="$1"
infilename=$(basename "$infile")
prg_name="${infilename%%.*}"
outdir="$2"
outfile="${outdir}/${prg_name}.deduplicated.fa"

# snippet taken from https://www.biostars.org/p/3003
STATS=$(fastx_collapser -i "$infile" -o "$outfile" -v 2>&1)

printf "====\n%s stats:\n%s\n" "$prg_name" "$STATS" >&2
