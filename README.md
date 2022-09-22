# [IN DEVELOPMENT]

# yEvo Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.14.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![Conda](https://anaconda.org/conda-forge/terraform-provider-github/badges/installer/conda.svg)](https://docs.conda.io/en/latest/miniconda.html)

Snakemake pipeline for analyzing sequencing data from the [yEvo](https://yevo.org/) project.

## Installation

1. Make sure you have [conda](https://docs.conda.io/en/latest/miniconda.html) installed. 

2. Install Mamba to facilitate snakemake installation, as recommended in the [Snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba).

```
$ conda install -n base -c conda-forge mamba
```

3. Clone this repo, then create and activate the provided [environment](./environment.yml):

```
$ git clone https://github.com/dunhamlab/yevo_pipeline.git \
    && cd yevo_pipeline/ \
    && mamba env create -f environment.yml \
    && conda activate yevo_pipeline_env
```
