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

3. Clone this repo:

```
$ git clone https://github.com/dunhamlab/yevo_pipeline.git
```

4. Create the provided [environment](./environment.yml):

```
$ cd yevo_pipeline/ && mamba env create -f environment.yml
```

5. Activate the new conda environment:

```
$ conda activate yevo_pipeline_env
```

6. Downloaded required pipeline inputs and test sequencing data:

```
$ ./scripts/download_test_data.sh
```

7. Generate the `run_pipeline.sh` script using the included utility script:

```
$ ./scripts/gen_run_script.sh
```

You're ready to run the pipeline!

## Running the Pipeline


After following the above installation instructions, run the pipeline on the provided test input files:

```
$ ./run_pipeline.sh
```

**NOTE:** be sure that you are in the repo's base directory with the `yevo_pipeline_env` conda environment activated.

To run this pipeline on your own sequencing data:

* in `run_pipeline.sh`, modify `FASTQ_DIR` to point to your raw fastq data dir
* in `run_pipeline.sh`, modify `OUTPUT_DIR` to point to your desired output dir

The paths to reference genome, ancestor, and annotation files can be configured using the [config/config.yml](config/config.yml) file.



