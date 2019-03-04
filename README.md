## Prerequisites
You need to have Singularity and snakemake available in your path.

A `Pipfile` is provided that can be used to create a python virtual environment
with `snakemake` in it. To activate, from the project root directory run

```sh
pipenv shell
```

You also need to have a singularity image available that can run the ONT
basecaller Guppy. I cannot make an image for this available as it is not open-source
software.

## Setup

### Data

This pipeline requires that your data be in a particular location with a specified
naming system for the directory names.  

Firstly, create a directory within the project root called `data`.  

Within `data` you should create a directory for each "region"/"project". Then,
within each region/project there should be a directory for each nanopore run.  

In each nanopore run directory there should be a directory called `f5s` and all
of the fast5 files for that run should be somewhere under this directory.  

As an example, if I had a region/project called "madagascar" and 2 runs called "a"
and "b", then my `data` directory would look like this

```
├── data
│   └── madagascar
│       ├── a
│       │   └── f5s
│       └── b
│           └── f5s
```

### Sample sheet
Another required part of the pipeline is a sample sheet. The location and name
for this sheet should be specified under the key [`samples`](https://github.com/mbhall88/head_to_head_pipeline/blob/1118b0b6ffc158b1de5792651775d64f5f0d7562/config.yaml#L1),
but it must be a CSV file. The default for this sample sheet is in `docs/samples.csv`.  

The sample sheet structure will be validated when running the pipeline, so you
need to ensure it is formatted correctly or the pipeline wont run.  

The sample sheet must have 6 columns with the following information in each row:
1. Region - the name of the region for that sample. This must also be the name
of the region directory in the `data` directory as [discussed above](#data) (`madagascar` would be an example).
2. Nanopore run ID - the name for the nanopore run ID this sample comes from. In the
example from [above](#data) this would be `a` or `b`.
3. Sample ID - this must be a unique identifier for this sample.
4. GUUID - a unique identifier for the Illumina data for this sample (not required).
5. Nanopore barcode - the nanopore barcode number used for this sample (if multiplexed).
This must be entered as barcodeXX where XX is the number i.e 01, 04 etc. So for
a sample that was barcoded with number 4 you would enter `barcode04`. If the
sample was not multiplexed, then leave this blank.
6. Flowcell version - the chemistry version of the flowcell. This defaults to
R9.4.1. Currently, only R9.4 and R9.5 are supported.

## Run
*NOTE: at this stage these instructions are for running on an LSF cluster system*  

To execute the pipeline, make sure snakemake and Singularity are available in your
path. Make sure you are in the pipeline root directory and run the following

```sh
bash analysis/scripts/submit_lsf.sh snakemake_master_process
```
