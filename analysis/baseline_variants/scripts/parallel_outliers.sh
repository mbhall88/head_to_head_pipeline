#!/usr/bin/env bash
compass_dir="/hps/nobackup/research/zi/projects/tech_wars/analysis/baseline_variants/illumina/gvcfs/"
mask="/hps/nobackup/research/zi/projects/tech_wars/analysis/baseline_variants/resources/compass-mask.bed"
indir="$1"
outdir="$2"
threads="$3"

for outlier in "mada_2-46" "R27006" "R13303" "mada_1-34" "mada_1-50"
do
echo "Running concordance for $outlier"
a="${compass_dir}/${outlier}.compass.vcf.gz"
outlier_dir="${outdir}/$outlier"
mkdir -p "$outlier_dir"
fd -j "$threads" -HI \
    -x python -f -a "$a" -b '{}' -m "$mask" -c "${outlier_dir}/"'{/.}'".concordance.csv" -j "${outlier_dir}/"'{/.}'".concordance.json" \; \
    'snps.filtered.bcf$' "$indir"
done
