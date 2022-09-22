
#
# yEvo Pipeline
#

import os
from datetime import datetime

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Define Constants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# discover input files
SAMPLES, READS = glob_wildcards(f"{config['fastq_dir']}/{{sample}}_{{read}}_001.fastq.gz")
SAMPLES = list(set(SAMPLES))
READS = list(set(READS))

# create a new timestamped output directory for every pipeline run
OUTPUT_DIR = f"results/{datetime.now().strftime('%Y%m%d_%H%M%S')}"

# Project name and date for bam header
SEQID='yEvo_hackathon_align'


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Begin Pipeline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# https://snakemake.readthedocs.io/en/v7.14.0/tutorial/basics.html#step-7-adding-a-target-rule 
rule all:
    input:
        f'{OUTPUT_DIR}/DONE.txt'


# ~~~~~~~~~~~~~~~~~~~~~~~~~~ Set Up Reference Files ~~~~~~~~~~~~~~~~~~~~~~~~~ #

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

        
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Run QC on Raw Reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

rule run_fastqc_all:
    input:
        config['fastq_dir']
    output:
        f'{OUTPUT_DIR}/02_fastqc/all_samples/stdin_fastqc.html',
        f'{OUTPUT_DIR}/02_fastqc/all_samples/stdin_fastqc.zip',
    conda:
        'envs/main.yml'
    shell:
        f'zcat {input}/*.fastq.gz | fastqc stdin --outdir={OUTPUT_DIR}/02_fastqc/all_samples/'
        

rule run_fastqc_persample:
    input:
        f"{config['fastq_dir']}/{{sample}}_R1_001.fastq.gz",
        f"{config['fastq_dir']}/{{sample}}_R2_001.fastq.gz",
    output:
        f"{OUTPUT_DIR}/02_fastqc/{{sample}}/stdin_fastqc.html",
        f"{OUTPUT_DIR}/02_fastqc/{{sample}}/stdin_fastqc.zip"
    conda:
        'envs/main.yml'
    shell:
        f"zcat {{input}} | fastqc stdin --outdir={OUTPUT_DIR}/02_fastqc/{{wildcards.sample}}"


# ~~~~~~~~~~~~~~~~~~~~~~~~ Perform Initial Alignment ~~~~~~~~~~~~~~~~~~~~~~~~ #

rule align_reads:
    input:
        rules.create_bwa_index.output,
        r1=f"{config['fastq_dir']}/{{sample}}_R1_001.fastq.gz",
        r2=f"{config['fastq_dir']}/{{sample}}_R2_001.fastq.gz",
    output:
        f"{OUTPUT_DIR}/03_init_alignment/{{sample}}_R1R2.sam"
    conda:
        'envs/main.yml'
    shell:
        f"bwa mem -R '@RG\tID:{SEQID}\tSM:{{wildcards.sample}}\tLB:1' {{rules.copy_fasta.output}} {{input.r1}} {{input.r2}} > {{output}}"
        

rule samtools_view:
    input:
        rules.align_reads.output
    output:
        f"{OUTPUT_DIR}/03_init_alignment/{{sample}}_R1R2.bam"
    conda:
        'envs/main.yml'
    shell:
        "samtools view -bS {input} -o {output}"


rule samtools_sort_one:
    input:
        rules.samtools_view.output
    output:
        f"{OUTPUT_DIR}/03_init_alignment/{{sample}}_R1R2_sort.bam"
    conda:
        'envs/main.yml'
    shell:
        "samtools sort {input} -o {output}"
        

rule samtools_index_one:
    input:
        rules.samtools_sort_one.output
    output:
        f"{OUTPUT_DIR}/03_init_alignment/{{sample}}_R1R2_sort.bam.bai"
    conda:
        'envs/main.yml'
    shell:
        "samtools index {input}"
        

rule samtools_flagstat:
    input:
        rules.samtools_sort_one.output
    output:
        f"{OUTPUT_DIR}/03_init_alignment/{{sample}}_R1R2_sort_flagstat.txt"
    conda:
        'envs/main.yml'
    shell:
        "samtools flagstat {input} > {output}"


rule picard_mark_dupes:
    input:
        rules.samtools_sort_one.output
    output:
        bam=f"{OUTPUT_DIR}/04_picard/{{sample}}_R1R2_comb.MD.bam",
        metrics=f"{OUTPUT_DIR}/04_picard/{{sample}}_comb_R1R2.sort_dup_metrix"
    conda:
        'envs/main.yml'
    shell:
        "picard MarkDuplicates --INPUT {input} --OUTPUT {output.bam} --METRICS_FILE {output.metrics} --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY LENIENT"


rule picard_read_groups:
    input:
        rules.picard_mark_dupes.output.bam
    output:
        f"{OUTPUT_DIR}/04_picard/{{sample}}_R1R2_comb.RG.MD.bam",
    conda:
        'envs/main.yml'
    shell:
        f"picard AddOrReplaceReadGroups --INPUT {{input}} --OUTPUT {{output}} --RGID {SEQID} --RGLB 1 --RGPU 1 --RGPL illumina --RGSM {{wildcards.sample}} --VALIDATION_STRINGENCY LENIENT"







rule finish:
    input:
        rules.run_fastqc_all.output,
        expand(rules.run_fastqc_persample.output, sample=SAMPLES),
        expand(rules.picard_read_groups.output, sample=SAMPLES),
        expand(rules.samtools_index_one.output, sample=SAMPLES),
        expand(rules.samtools_flagstat.output, sample=SAMPLES),
    output:
        f'{OUTPUT_DIR}/DONE.txt'
    shell:
        'touch {output}'
