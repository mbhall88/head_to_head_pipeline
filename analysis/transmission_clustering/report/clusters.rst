An interactive plot of the {{ snakemake.wildcards.caller }} clusters produced by samples within a certain distance of
each other. To select the distance threshold, a linear regression is fitted to the pairwise relationship of Illumina
(COMPASS) and Nanopore ({{ snakemake.wildcards.caller }}) SNP distances. We then set the regression coefficient (*x*)
to the Illumina SNP threshold of {{ snakemake.wildcards.threshold }} in this relationship to obtain the response variable
(*y*; Nanopore distance threshold). Importantly, we only fit this model based on pairs of samples within a distance of
{{ snakemake.config["adaptive_threshold"] }} SNPs of each other. See the log file for the Nanopore threshold used.

The first plot shows the true (Illumina) clusters with the nodes coloured according to their precision (outer colour)
and recall (inner colour). For an explanation of how these metrics are defined and calculated, see `this github comment`_.

The second plot shows the same as the second, but the nodes are coloured by their cluster's average precision and recall
values.

.. _`this github comment`: https://github.com/mbhall88/head_to_head_pipeline/issues/65#issuecomment-797895910