The relationship of the distance between all pairs of samples based on COMPASS VCF calls (X-axis) and mixed COMPASS-``bcftools`` calls (Y-axis).
The black, dashed line indicates the relationship we would expect if the distance between a pair of samples were the same for both approaches.
The blue line indicates the line of best fit based on fitting a robust linear regression model to the data. The inset gives a closer look at the
relationship for all sample pairs where the COMPASS distance is less than or equal to {{ snakemake.params.inset_threshold }} SNPs. The legend
indicates the linear equations for the lines.
*Note: to prevent model skew, we do not include self-distance pairs.*