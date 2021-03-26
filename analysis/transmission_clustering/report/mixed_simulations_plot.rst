This plot shows the results of simulation varying ratios (X-axis) of Nanopore/Illumina mixtures. For example, a ratio of 0.75
means that, for clustering, we use Nanopore SNP calls for 80% of the samples, and Illumina for the other 20%. The
different thresholds (subplots) indicate the cutoff for defining two samples as "connected". The Y-axis depicts the
Sample-Averaged Cluster Recall and Precision (SACR/SACP) and Cluster Number Ratio (CNR) distributions for {{ snakemake.params.num_simulations }}
simulations - i.e. the number of times we randomly split the samples for the given ratio.

SACR: We define recall for each sample in the "truth" (Illumina/COMPASS) clustering, to be the proportion of the members
of the sample's truth cluster that appear in the mixed technology clustering. We average all of the values to give the
SACR. This metric is used to illustrate whether we miss samples from clustering. A value of 1.0 means that all samples
that are supposed to be clustered together *are* clustered together in the mixed technology clustering.

SACP: Similar to the SACR, we define precision to be the proportion of members in the sample's mixed technology clustering
that also appear in its truth cluster. This metric is designed to illustrate when we add samples to clusters that shouldn't
be in that cluster. A value of 1.0 means mixed technology clustering adds no excess samples.

CNR: The number of clusters in the COMPASS/Illumina clustering, divided by the number of clusters in the mixed technology
clustering. This metric is intended to give an idea of how many false positive clusters the mixed technology identifies.