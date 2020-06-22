#!/usr/bin/env bash
set -u
infile="$1"
outfile="$2"
logfile="$3"
linsi --thread -1 "$infile" > "$outfile" 2>> "$logfile"