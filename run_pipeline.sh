#!/bin/bash

# activate conda env
source activate yevo_pipeline_env

# run the pipeline
snakemake --cores --snakefile "$SNAKE_FILE" --configfile "$CONFIG_FILE" \
--config output_dir=$OUTPUT_DIR fastq_dir=$FASTQ_DIR anc_dir=$ANC_DIR ref_fasta=$REF_FASTA\
--use-conda --conda-prefix="$HOME/.snakemake/conda"

# success
echo -e "\nDONE!\n"
