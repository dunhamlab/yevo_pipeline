#!/bin/bash

#
# configure file paths
#

# absolute path to the cloned yevo_pipeline repo
REPO_DIR="{{BASE_DIR}}"

# absolute path to fastq.gz location (e.g. <sample>_R1_001.fastq.gz files)
FASTQ_DIR="$REPO_DIR/data/fastq"

# absolute path to desired pipeline output directory (will be created by snakemake)
OUTPUT_DIR="$REPO_DIR/results/test_run"

#
# run pipeline
#

# go to working directory
cd $REPO_DIR

# construct pipeline file paths
LOG_FILE="$OUTPUT_DIR/00_logs/yevo_pipeline.log"
SNAKE_FILE="$REPO_DIR/workflow/Snakefile.py"
CONFIG_FILE="$REPO_DIR/config/config.yml"

# activate conda env
source activate yevo_pipeline_env

# run the pipeline
snakemake --cores --snakefile $SNAKE_FILE --configfile $CONFIG_FILE \
    --config output_dir=$OUTPUT_DIR fastq_dir=$FASTQ_DIR \
    --use-conda --conda-prefix="$HOME/.snakemake/conda"

# success
echo -e "\nDONE!\n"
