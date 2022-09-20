
#
# yEvo Pipeline
#

import os
from datetime import datetime

# create a new timestamped output directory for every pipeline run
OUTPUT_DIR = f"results/{datetime.now().strftime('%Y%m%d_%H%M%S')}"

# https://snakemake.readthedocs.io/en/v7.14.0/tutorial/basics.html#step-7-adding-a-target-rule 
rule all:
    input:
        f'{OUTPUT_DIR}/DONE.txt'




#
# copy the supplied reference genome fasta to the pipeline output directory for reference
#
rule copy_fasta:
    input:
        config['ref_fasta']
    output:
        f"{OUTPUT_DIR}/01_ref_files/{os.path.basename(config['ref_fasta'])}"
    shell:
        "cp {input} {output}"

#
# create a BWA index from the copied fasta reference genome
#
rule create_bwa_index:
    input:
        rules.copy_fasta.output
    output:
        f"{rules.copy_fasta.output}.amb",
        f"{rules.copy_fasta.output}.ann",
        f"{rules.copy_fasta.output}.bwt",
        f"{rules.copy_fasta.output}.pac",
        f"{rules.copy_fasta.output}.sa",
    conda:
        'envs/main.yml'
    shell:
        "bwa index {input}"

        

#
# run fastqc to generate quality control reports on the raw reads
#
rule run_fastqc:
    input:
        config['fastq_dir']
    output:
        f'{OUTPUT_DIR}/02_fastqc/stdin_fastqc.html',
        f'{OUTPUT_DIR}/02_fastqc/stdin_fastqc.zip',
    conda:
        'envs/main.yml'
    shell:
        f'zcat {input}/*.fastq.gz | fastqc stdin --outdir={OUTPUT_DIR}/02_fastqc'
        

        
        
        
        

rule finish:
    input:
        rules.run_fastqc.output,
        rules.create_bwa_index.output,
    output:
        f'{OUTPUT_DIR}/DONE.txt'
    conda:
        'envs/main.yml'
    shell:
        'touch {output}'
