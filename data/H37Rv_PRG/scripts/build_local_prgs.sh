#!/usr/bin/env sh
set -u
msa="$1"
outdir="$2"
out_prefix="${outdir}/${$(basename "$infile")%.*}"
nesting_lvl="$3"
match_len="$4"
logfile="$5"

make_prg prg_from_msa --max_nesting "$nesting_lvl" \
    --min_match_length "$match_len" \
    --prefix "$out_prefix" \
    "$msa" 2>> "$logfile"
cat "$out_prefix".log >> "$logfile"
rm "$out_prefix".log "$out_prefix".gfa