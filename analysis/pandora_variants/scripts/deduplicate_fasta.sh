#!/usr/bin/env bash
set -eu

infile="$1"
infilename=$(basename "$infile")
outdir="$2"
outfile="${outdir}/${infilename%%.*}.deduplicated.fa"

# snippet taken from https://www.biostars.org/p/3003/#3008
sed -e '/^>/s/$/@/' -e 's/^>/#/' "$infile"  |\
tr -d '\n' | tr "#" "\n" | tr "@" "\t" |\
sort -u -t '  ' -f -k 2,2  |\
sed -e 's/^/>/' -e 's/\t/\n/' > "$outfile"
