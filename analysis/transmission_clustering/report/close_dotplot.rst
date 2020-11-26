An interactive dotplot for a subset of distance matrices produced by |snp-dists|_ for
{{ snakemake.params.xname }} and {{ snakemake.params.ynames }} where the
{{ snakemake.params.xname }} SNP distance is no greater than {{ snakemake.params.threshold }}.
Each point is the SNP distance between a pair of samples. The x-axis is the distance
between pairs based on SNP calls from {{ snakemake.params.xname }} (Illumina), and the y-axis is
for {{ snakemake.params.ynames }} (Nanopore) calls.

.. |snp-dists| replace:: ``snp-dists``
.. _snp-dists: https://github.com/tseemann/snp-dists