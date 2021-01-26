#!/usr/bin/env bash
vcf="$1"
matrix="$2"

printf "CHROM,POS," >"$matrix"
(bcftools query -l "$vcf" | tr '\n' ',' | sed 's/,$/\n/g') >>"$matrix"
(bcftools query -f '%CHROM,%POS,[%GT,]\n' "$vcf" | sed 's/,$//g') >>"$matrix"
