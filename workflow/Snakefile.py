
#
# yEvo Pipeline
#

import os
from datetime import datetime

# create a new timestamped output directory for every pipeline run
OUTPUT_DIR = f"results/{datetime.now().strftime('%Y%m%d_%H%M%S')}"

# discover input files
SAMPLES, READS = glob_wildcards(f"{config['fastq_dir']}/{{sample}}_{{read}}_001.fastq.gz")
SAMPLES = list(set(SAMPLES))
READS = list(set(READS))

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
        


#
# run fastqc to generate quality control reports on the raw reads
#
rule align_reads:
    input:
        rules.create_bwa_index.output,
        r1=f"{config['fastq_dir']}/{{sample}}_R1_001.fastq.gz",
        r2=f"{config['fastq_dir']}/{{sample}}_R2_001.fastq.gz",
    output:
        f"{OUTPUT_DIR}/03_bwa/{{sample}}.sam"
    conda:
        'envs/main.yml'
    shell:
        "bwa mem -R '@RG\tID:YEVO_SEQID\tSM:{wildcards.sample}\tLB:1' {rules.copy_fasta.output} {input.r1} {input.r2} > {output}"
       
        
        
 


        
        

rule finish:
    input:
        expand(rules.align_reads.output, sample=SAMPLES),
        rules.run_fastqc.output,
    output:
        f'{OUTPUT_DIR}/DONE.txt'
    shell:
        'touch {output}'
