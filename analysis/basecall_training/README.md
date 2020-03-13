Scripts and workflows for training a *M. tuberculosis* species-specific nanopore 
basecalling model.

## Prerequisites
To run the whole pipeline you will need the following programs installed:
-   [Singularity (v3)][singularity]
-   [Conda][conda]
-   [Snakemake][snakemake]
-   [ripgrep][ripgrep]

## Setup
If you cloned this repository with `git`, without also cloning the submodules, you will 
additionally need to fetch the files for the submodules this workflow relies on.

```bash
git submodule update --init --recursive
```

[singularity]: https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps
[conda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/
[snakemake]: https://snakemake.readthedocs.io/en/stable/
[ripgrep]: https://github.com/BurntSushi/ripgrep
