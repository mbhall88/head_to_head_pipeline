An interactive plot of the {{ snakemake.wildcards.caller }} clusters produced by samples within a certain distance of
each other. To select the distance threshold, a linear regression is fitted to the pairwise relationship of Illumina
(COMPASS) and Nanopore ({{ snakemake.wildcards.caller }}) SNP distances. We then set the regression coefficient (*x*)
to the Illumina SNP threshold of {{ snakemake.wildcards.threshold }} in this relationship to obtain the response variable
(*y*; Nanopore distance threshold). Importantly, we only fit this model based on pairs of samples within a distance of
{{ snakemake.config["adaptive_threshold"] }} SNPs of each other. See the log file for the Nanopore threshold used.