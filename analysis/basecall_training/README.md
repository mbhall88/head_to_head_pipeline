Scripts and workflows for training a *M. tuberculosis* species-specific nanopore
basecalling model.

[TOC]:#

# Table of Contents
- [Prerequisites](#prerequisites)
- [Setup](#setup)
- [Usage](#usage)


## Prerequisites

To run the whole pipeline you will need the following programs installed:
- [Singularity (v3)][singularity]
- [Conda][conda]
- [Snakemake][snakemake]
- [ripgrep][ripgrep]

## Setup

If you cloned this repository with `git`, without also cloning the submodules, you will
additionally need to fetch the files for the submodules this workflow relies on.

```bash
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

```bash
# where the CUDA binaries and libraries are stored on my cluster
cuda_bins="/usr/local/cuda/bin"
cuda_libs="/usr/local/cuda/lib64/"
# on my cluster I need to manually mount the CUDA directories into the container
bash scripts/submit_lsf.sh --singularity-args "--nv --bind $cuda_bins --bind $cuda_libs"
```


[singularity]: https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps
[conda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/
[snakemake]: https://snakemake.readthedocs.io/en/stable/
[ripgrep]: https://github.com/BurntSushi/ripgrep
[lsf-profile]: https://github.com/Snakemake-Profiles/snakemake-lsf

