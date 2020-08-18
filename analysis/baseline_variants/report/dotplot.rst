An interactive dotplot for the distance matrices produced by |snp-dists|_ for {{ snakemake.params.xname }} and {{ snakemake.params.yname }}.
Each point is the SNP distance between a pair of samples. The x-axis is the distance
between pairs based on SNP calls from {{ snakemake.params.xname }}, and the y-axis is
for {{ snakemake.params.yname }} calls.

.. |snp-dists| replace:: ``snp-dists``
.. _snp-dists: https://github.com/tseemann/snp-dists