<div align="center">
    <a href="#readme"><img src="./docs/img/yevo_pipeline_logo.png" width="300"></a>
</div>


# yEvo Pipeline

[![Snakemake](./docs/img/snakemake.svg)](https://snakemake.readthedocs.io)
[![Conda](./docs/img/conda.svg)](https://docs.conda.io/en/latest/miniconda.html)

Snakemake pipeline for analyzing sequencing data from the [yEvo](https://yevo.org/) project.

![yEvo Pipeline](../../raw/main/docs/img/pipeline.svg)

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

4. Downloaded required pipeline inputs and test sequencing data:

```
$ ./scripts/download_test_data.sh
```

## Running the Pipeline

TODO - coming soon!
