#!/usr/bin/env sh
set -eux
msa="$1"
out_prefix="$2"
nesting_lvl="$3"
match_len="$4"
logfile="$5"

make_prg prg_from_msa --max_nesting "$nesting_lvl" \
    --min_match_length "$match_len" \
    --prefix "$out_prefix" \
    "$msa" 2>> "$logfile"

full_outprefix="${out_prefix}.max_nest${nesting_lvl}.min_match${match_len}"
cat "${full_outprefix}.log" >> "$logfile"
rm -f "$full_outprefix".log "$full_outprefix".gfa "$full_outprefix".bin
