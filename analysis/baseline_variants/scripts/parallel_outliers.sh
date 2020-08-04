#!/usr/bin/env bash
compass_dir="/hps/nobackup/research/zi/projects/tech_wars/analysis/baseline_variants/illumina/gvcfs/"
mask="/hps/nobackup/research/zi/projects/tech_wars/analysis/baseline_variants/resources/compass-mask.bed"
pyscript="$1"
indir="$2"
outdir="$3"
threads="$4"

for outlier in "mada_2-46" "R27006" "R13303" "mada_1-34" "mada_1-50"
do
echo "Running concordance for $outlier"
b="${indir}/*/${outlier}.snps.filtered.bcf"
outlier_dir="${outdir}/$outlier"
mkdir -p "$outlier_dir"
fd -j "$threads" -HI \
    -x python "$pyscript" -f -a '{}' -b "$b" -m "$mask" -c "${outlier_dir}/"'{/.}'".concordance.csv" -j "${outlier_dir}/"'{/.}'".concordance.json" \; \
    'vcf.gz$' "$compass_dir"
done
