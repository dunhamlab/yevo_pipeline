#!/bin/bash

#$ -l mfree=8778M
#
# configure file paths
#

# absolute path to the cloned yevo_pipeline repo
REPO_DIR="/net/dunham/vol1/home/dennig2/yevo_pipeline"

# absolute path to fastq.gz location (e.g. <sample>_R1_001.fastq.gz files)
FASTQ_DIR="$REPO_DIR/data/fastq"

# absolute path to desired pipeline output directory (will be created by snakemake)
OUTPUT_DIR="$REPO_DIR/data/output"

#
# run pipeline
#

# construct pipeline file paths
LOG_FILE="$OUTPUT_DIR/00_logs/yevo_pipeline.log"
SNAKE_FILE="$REPO_DIR/workflow/Snakefile.py"
CONFIG_FILE="$REPO_DIR/config/config.yml"
SAMPLE_LIST="$OUTPUT_DIR/00_logs/00_sample_list.txt"

# go to working directory
cd $REPO_DIR

# activate conda env
source activate yevo_pipeline

# run the pipeline
snakemake --use-conda --conda-prefix="$HOME/.snakemake/conda" --cores 4 --snakefile $SNAKE_FILE --configfile $CONFIG_FILE --config output_dir=$OUTPUT_DIR fastq_dir=$FASTQ_DIR --rerun-incomplete --keep-going

# success
echo -e "\nDONE!\n"
