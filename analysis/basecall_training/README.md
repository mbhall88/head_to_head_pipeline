Scripts and workflows for training a *M. tuberculosis* species-specific nanopore
basecalling model.

[TOC]: #

# Table of Contents
- [Prerequisites](#prerequisites)
- [Setup](#setup)
- [Usage](#usage)
- [Method](#method)
- [Results](#results)


## Prerequisites

To run the whole pipeline you will need the following programs installed:
- [Singularity (v3)][singularity]
- [Conda][conda]
- [Snakemake][snakemake]
- [ripgrep][ripgrep]

## Setup

Recursively clone the repository

```sh
# git version >=2.13
git clone --recursive-submodules https://github.com/mbhall88/head_to_head_pipeline.git
# git version <=2.13
git clone --recursive https://github.com/mbhall88/head_to_head_pipeline.git

cd head_to_head_pipeline/analysis/basecall_training 
```

If you cloned this repository with `git`, without also cloning the submodules, you will
additionally need to fetch the files for the submodules this workflow relies on.

```sh
git submodule update --init --recursive
```

## Usage

When submitting on an LSF cluster you may wish to use the provided submission script
[`scripts/submit_lsf.sh`](scripts/submit_lsf.sh). This script assumes you are using the
[snakemake LSF profile][lsf-profile].

Additionally, the provided `lsf.yaml` configuration file sets out the LSF parameters to
use for submitting the model training rule to a multi-GPU cluster. Note that the host
(`m`) option may differ on your cluster. If using GPUs ensure you pass
`--singularity-args '--nv'` to `snakemake` to tell Singularity to use the local GPU
resources.

An example of how I submit this pipeline on my LSF cluster is

```sh
bash scripts/submit_lsf.sh --singularity-args "--nv"
```

Depending on the size of the data, this pipeline may take a week or two to run.

## Method

An overview of the workflow is below

![workflow overview](resources/rulegraph.png)

The bulk of the training-related rules are explained in further detail in the [taiyaki
walkthrough][walkthrough].

The additional steps for evaluating the pipeline involve mapping the reads that were set
aside for evaluation to the truth assemblies and then calculating the accuracy metrices
for this. See the report notebook (below) for more information about these steps.

Versions for software used can be found in [`config.yaml`](config.yaml).

## Results

All of the code, with inline plots, can be found at [`report/processed-report.ipynb`](report/processed-report.ipynb). Or alternatively, in a nice rendered format [here][nbviewer].



[singularity]: https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps
[conda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/
[snakemake]: https://snakemake.readthedocs.io/en/stable/
[ripgrep]: https://github.com/BurntSushi/ripgrep
[lsf-profile]: https://github.com/Snakemake-Profiles/snakemake-lsf
[walkthrough]: https://github.com/nanoporetech/taiyaki/blob/master/docs/walkthrough.rst
[nbviewer]: https://nbviewer.jupyter.org/github/mbhall88/head_to_head_pipeline/blob/master/analysis/basecall_training/report/processed-report.ipynb
[blast]: https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity#blast-identity
