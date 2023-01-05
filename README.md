### Paper

> Hall, M. B. et al. Evaluation of Nanopore sequencing for Mycobacterium tuberculosis drug susceptibility testing and outbreak investigation: a genomic analysis. *The Lancet Microbe* 0, (2022) doi: [10.1016/S2666-5247(22)00301-9][doi].

[doi]: https://doi.org/10.1016/S2666-5247(22)00301-9

---

This repository holds the pipelines/scripts used for our paper analysing Illumina and
Nanopore for *M.tuberculosis* drug resistance calling and transmission clustering.

For people wanting to analyse their Nanopore data in the same manner as we did in this paper, we would suggest using https://github.com/mbhall88/tbpore, which is a python program that runs the drug resistance prediction and clustering (with a smaller decontamination database) components of this pipeline. It is actively maintained and much easier to use. 

All pipelines require the following dependencies to be installed:
- [Snakemake](https://snakemake.github.io/)
- [Conda](https://docs.conda.io/en/latest/) (and
  [Mamba](https://github.com/mamba-org/mamba))
- [Singularity](https://sylabs.io/docs)
- The Python library [`pandas`](https://pandas.pydata.org/)

See subdirectories for more specific information about different pipelines. They are
nested according to their dependence on the outputs of each pipeline.

- [Quality Control](data/QC)
  - [Assembly](analysis/assembly)
  - [Baseline variant analysis](analysis/baseline_variants)
    - [Transmission clustering](analysis/transmission_clustering)
  - [Drug Resistance Prediction](analysis/resistance_prediction)

The following pipelines are not relevant to the work in the final paper.

- [H37Rv PRG construction](data/H37Rv_PRG)
- [Pandora variant analysis](analysis/pandora_variants)

# Data availability

All data is submitted under the Project accession **PRJEB49093**.

The accessions and all relevant sample metadata for this study can be found at <https://doi.org/10.6084/m9.figshare.19304648>.

The raw Nanopore data is available to download from: <https://ftp.ebi.ac.uk/pub/databases/ont_tb_eval2022/>. See the sample metadata file for mappings between samples and the relevant Nanopore runs and barcode numbers.
