This repository holds the pipelines/scripts used for our project analysing Illumina and
Nanopore for Mtb drug resistance calling and tranmission clustering.

It is currently in progress.

All pipelines require the following dependencies to be installed:
- [Snakemake](https://snakemake.github.io/)
- [Conda](https://docs.conda.io/en/latest/) (and
  [Mamba](https://github.com/mamba-org/mamba))
- [Singularity](https://sylabs.io/docs)
- The Python library [`pandas`](https://pandas.pydata.org/)

See subdirectories for more specific information about different pipelines.

- [Quality Control](data/QC)
- [Assembly](analysis/assembly)
- [Baseline variant analysis](analysis/baseline_variants)
- [Transmission clustering](analysis/transmission_clustering)
- [Drug Resistance Prediction](analysis/resistance_prediction)

The following pipelines are not relevant to the work in the final paper.

- [H37Rv PRG construction](data/H37Rv_PRG)
- [Pandora variant analysis](analysis/pandora_variants)

