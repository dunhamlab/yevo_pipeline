
#
# yEvo Pipeline
#

import os
import json
from datetime import datetime


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Define Constants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# discover input files using path from run config
SAMPLES = list(set(glob_wildcards(f"{config['fastq_dir']}/{{sample}}_R1_001.fastq.gz").sample))

# read output dir path from run config
OUTPUT_DIR = config['output_dir']

# Project name and date for bam header
SEQID='yevo_pipeline_align'

REF_DIR = config['ref_fasta']


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Begin Pipeline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# https://snakemake.readthedocs.io/en/v7.14.0/tutorial/basics.html#step-7-adding-a-target-rule 
rule all:
    input:
        f'{OUTPUT_DIR}/DONE.txt'


# ~~~~~~~~~~~~~~~~~~~~~~~~~~ Set Up Reference Files ~~~~~~~~~~~~~~~~~~~~~~~~~ #

#
# export the current run configuration in JSON format
#
rule export_run_config:
    output:
        path=f"{OUTPUT_DIR}/00_logs/00_run_config.json"
    run:
        with open(output.path, 'w') as outfile:
            json.dump(dict(config), outfile, indent=4)


#
# make a list of discovered samples
#
rule list_samples:
    output:
        f"{OUTPUT_DIR}/00_logs/00_sample_list.txt"
    shell:
        "echo -e '{}' > {{output}}".format('\n'.join(SAMPLES))


#
# copy the supplied reference genome fasta to the pipeline output directory for reference
#
rule copy_fasta:
    input:
        config['ref_fasta']
    output:
        f"{OUTPUT_DIR}/01_ref_files/reference"
    shell:
        "cp {input} {output} && chmod 777 {output}"

rule copy_fasta_file:
    input:
        ref=f"{REF_DIR}"
    output:
        f"{OUTPUT_DIR}/01_ref_files/reference.fasta"   
    shell:
        "cp {input.ref} {output} && chmod 777 {output}"

rule index_fasta_file:
    input:
        rules.copy_fasta_file.output
    output:
        f"{rules.copy_fasta_file.output}.fai"
    conda:
        'envs/main.yml'
    shell:
        "samtools faidx {input}"

rule create_ref_dict:
    input:
        rules.copy_fasta_file.output
    output:
        f"{OUTPUT_DIR}/01_ref_files/reference.dict"     
    conda:
        'envs/main.yml'
    shell:
        "picard CreateSequenceDictionary -R {input} -O {output} && chmod 777 {output}"

# create a BWA index from the copied fasta reference genome
#
rule create_bwa_index:
    input:
        rules.copy_fasta_file.output
    output:
        f"{rules.copy_fasta_file.output}.amb",
        f"{rules.copy_fasta_file.output}.ann",
        f"{rules.copy_fasta_file.output}.bwt",
        f"{rules.copy_fasta_file.output}.pac",
        f"{rules.copy_fasta_file.output}.sa",
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
        f'cat {input}/*.fastq.gz | zcat | fastqc stdin --outdir={OUTPUT_DIR}/02_fastqc/all_samples/'
        

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
        f"cat {{input}} | zcat | fastqc stdin --outdir={OUTPUT_DIR}/02_fastqc/{{wildcards.sample}}"


# ~~~~~~~~~~~~~~~~~~~~~~~~ Perform Initial Alignment ~~~~~~~~~~~~~~~~~~~~~~~~ #

rule align_reads:
    input:
        rules.create_bwa_index.output,
        r1=f"{config['fastq_dir']}/{{sample}}_R1_001.fastq.gz",
        r2=f"{config['fastq_dir']}/{{sample}}_R2_001.fastq.gz",
    output:
        f"{OUTPUT_DIR}/03_init_alignment/{{sample}}/{{sample}}_R1R2_sort.bam"
    conda:
        'envs/main.yml'
    shell:
        r"""bwa mem -R '@RG\tID:""" + SEQID + r"""\tSM:""" + '{wildcards.sample}' + r"""\tLB:1'""" + ' {rules.copy_fasta_file.output} {input.r1} {input.r2} | samtools sort -o {output} -'


rule samtools_index_one:
    input:
        rules.align_reads.output
    output:
        f"{OUTPUT_DIR}/03_init_alignment/{{sample}}/{{sample}}_R1R2_sort.bam.bai"
    conda:
        'envs/main.yml'
    shell:
        "samtools index {input}"
        

rule samtools_flagstat:
    input:
        bam=rules.align_reads.output,
        idx=rules.samtools_index_one.output
    output:
        f"{OUTPUT_DIR}/03_init_alignment/{{sample}}/{{sample}}_R1R2_sort_flagstat.txt"
    conda:
        'envs/main.yml'
    shell:
        "samtools flagstat {input.bam} > {output}"


rule picard_mark_dupes:
    input:
        bam=rules.align_reads.output,
        idx=rules.samtools_index_one.output,
        dct=rules.create_ref_dict.output
    output:
        bam=f"{OUTPUT_DIR}/04_picard/{{sample}}/{{sample}}_comb_R1R2.MD.bam",
        metrics=f"{OUTPUT_DIR}/04_picard/{{sample}}/{{sample}}_comb_R1R2.sort_dup_metrix"
    conda:
        'envs/main.yml'
    shell:
        "picard MarkDuplicates --INPUT {input.bam} --OUTPUT {output.bam} --METRICS_FILE {output.metrics} --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY LENIENT"


rule picard_read_groups:
    input:
        rules.picard_mark_dupes.output.bam
    output:
        f"{OUTPUT_DIR}/04_picard/{{sample}}/{{sample}}_comb_R1R2.RG.MD.sort.bam",
    conda:
        'envs/main.yml'
    shell:
        f"picard AddOrReplaceReadGroups --INPUT {{input}} --OUTPUT /dev/stdout --RGID {SEQID} --RGLB 1 --RGPU 1 --RGPL illumina --RGSM {{wildcards.sample}} --VALIDATION_STRINGENCY LENIENT | samtools sort -o {{output}} -"


rule samtools_index_two:
    input:
        rules.picard_read_groups.output
    output:
        f"{OUTPUT_DIR}/04_picard/{{sample}}/{{sample}}_comb_R1R2.RG.MD.sort.bam.bai",
    conda:
        'envs/main.yml'
    shell:
        "samtools index {input}"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GATK Re-alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule gatk_register:
    input:
        f'workflow/envs/src/GenomeAnalysisTK-3.7-0-gcfedb67.tar.bz2'
    output:
        f'{OUTPUT_DIR}/05_gatk/gatk_3.7_registered.txt'
    conda:
        'envs/main.yml'
    shell:
        "gatk-register {input} > {output} && chmod 777 {output}"


rule gatk_realign_targets:
    input:
        fa=rules.copy_fasta_file.output,
        bam=rules.picard_read_groups.output,
        idx=rules.samtools_index_two.output,
        gatk=rules.gatk_register.output,
        faidx=rules.index_fasta_file.output
    output:
        f'{OUTPUT_DIR}/05_gatk/{{sample}}/{{sample}}_comb_R1R2.bam.intervals'
    conda:
        'envs/main.yml'
    shell:
        "GenomeAnalysisTK -T RealignerTargetCreator -R {input.fa} -I {input.bam} -o {output}"


rule gatk_realign_indels:
    input:
        fa=rules.copy_fasta_file.output,
        intervals=rules.gatk_realign_targets.output,
        bam=rules.picard_read_groups.output,
        idx=rules.samtools_index_two.output        
    output:
        bam=f'{OUTPUT_DIR}/05_gatk/{{sample}}/{{sample}}_comb_R1R2.RG.MD.realign.bam',
        bai=f'{OUTPUT_DIR}/05_gatk/{{sample}}/{{sample}}_comb_R1R2.RG.MD.realign.bai'
    conda:
        'envs/main.yml'
    shell:
        "GenomeAnalysisTK -T IndelRealigner -R {input.fa} -I {input.bam} -targetIntervals {input.intervals} -o {output.bam}"


rule samtools_sort_three:
    input:
        rules.gatk_realign_indels.output.bam
    output:
        f'{OUTPUT_DIR}/05_gatk/{{sample}}/{{sample}}_comb_R1R2.RG.MD.realign.sort.bam',
    conda:
        'envs/main.yml'
    shell:
        "samtools sort {input} -o {output}"


rule samtools_index_three:
    input:
        rules.samtools_sort_three.output
    output:
        f'{OUTPUT_DIR}/05_gatk/{{sample}}/{{sample}}_comb_R1R2.RG.MD.realign.sort.bam.bai',
    conda:
        'envs/main.yml'
    shell:
        "samtools index {input}"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~ Begin Variant Calling ~~~~~~~~~~~~~~~~~~~~~~~~~~ #

rule bcftools_pileup:
    input:
        bam=rules.samtools_sort_three.output,
        idx=rules.samtools_index_three.output
    output:
        f'{OUTPUT_DIR}/06_variant_calling/{{sample}}/{{sample}}_samtools_AB.vcf',
    conda:
        'envs/main.yml'
    shell:
        "bcftools mpileup --ignore-RG -Ou -ABf {rules.copy_fasta_file.output} {input.bam} | bcftools call -vmO v -o {output}"


rule freebayes:
    input:
        bam=rules.samtools_sort_three.output,
        idx=rules.samtools_index_three.output
    output:
        f'{OUTPUT_DIR}/06_variant_calling/{{sample}}/{{sample}}_freebayes_BCBio.vcf',
    conda:
        'envs/main.yml'
    shell:
        "freebayes -f {rules.copy_fasta_file.output} --pooled-discrete --pooled-continuous --report-genotype-likelihood-max --allele-balance-priors-off --min-alternate-fraction 0.1 {input.bam} > {output}"


rule lofreq:
    input:
        bam=rules.samtools_sort_three.output,
        idx=rules.samtools_index_three.output,
        ancidx=f'{OUTPUT_DIR}/05_gatk/anc/anc_comb_R1R2.RG.MD.realign.sort.bam.bai',
        ancbam=f'{OUTPUT_DIR}/05_gatk/anc/anc_comb_R1R2.RG.MD.realign.sort.bam'
    output:
        normal=f'{OUTPUT_DIR}/06_variant_calling/{{sample}}/{{sample}}_lofreq_normal_relaxed.vcf.gz',
        tumor=f'{OUTPUT_DIR}/06_variant_calling/{{sample}}/{{sample}}_lofreq_tumor_relaxed.vcf.gz',
        somatic=f'{OUTPUT_DIR}/06_variant_calling/{{sample}}/{{sample}}_lofreq_somatic_final.snvs.vcf.gz',
    conda:
        'envs/main.yml'
    shell:
        f"lofreq somatic -n {{input.ancbam}} -t {{input.bam}} -f {{rules.copy_fasta_file.output}} -o {OUTPUT_DIR}/06_variant_calling/{{wildcards.sample}}/{{wildcards.sample}}_lofreq_"


rule unzip_lofreq:
    input:
        normal=rules.lofreq.output.normal,
        tumor=rules.lofreq.output.tumor,
        somatic=rules.lofreq.output.somatic,
    output:
        normal=rules.lofreq.output.normal.replace('.vcf.gz', '.vcf'),
        tumor=rules.lofreq.output.tumor.replace('.vcf.gz', '.vcf'),
        somatic=rules.lofreq.output.somatic.replace('.vcf.gz', '.vcf'),
    conda:
        'envs/main.yml'
    shell:
        "bgzip -d {input.normal} && bgzip -d {input.tumor} && bgzip -d {input.somatic}"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~ Begin Filtering Steps ~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#
# filter samtools results by ancestor
#
rule anc_filter_samtools:
    input:
        sample=rules.bcftools_pileup.output,
        anc=f'{OUTPUT_DIR}/06_variant_calling/anc/anc_lofreq_normal_relaxed.vcf'
    output:
        f'{OUTPUT_DIR}/07_filtered/{{sample}}/{{sample}}_samtools_AB_AncFiltered.vcf',    
    conda:
        'envs/main.yml'
    shell:
        f"bedtools intersect -v -header -a {{input.sample}} -b {{input.anc}} > {{output}}"


#
# filter freebayes results by ancestor
#
rule anc_filter_freebayes:
    input:
        sample=rules.freebayes.output,
        anc=f'{OUTPUT_DIR}/06_variant_calling/anc/anc_freebayes_BCBio.vcf'
    output:
        f'{OUTPUT_DIR}/07_filtered/{{sample}}/{{sample}}_freebayes_BCBio_AncFiltered.vcf',    
    conda:
        'envs/main.yml'
    shell:
        f"bedtools intersect -v -header -a {{input.sample}} -b {{input.anc}} > {{output}}"



#
# filter lofreq results by ancestor
#
rule anc_filter_lofreq:
    input:
        normal=rules.unzip_lofreq.output.normal,
        tumor=rules.unzip_lofreq.output.tumor,
    output:
        f'{OUTPUT_DIR}/07_filtered/{{sample}}/{{sample}}_lofreq_tumor_relaxed_AncFiltered.vcf',
    conda:
        'envs/main.yml'
    shell:
        f"bedtools intersect -v -header -a {{input.tumor}} -b {{input.normal}} > {{output}}"


# Many filtering steps
# MQ or MQM = Mapping quality
# QUAL = Metric that is specific for the variant caller that denotes confidence
# DP = Read depth
# DP[x] or SAF,SAR,SRF,SRR = Array with read depth for Fwd, Rev, and from which strand
# Filters by quality, mapping quality, read depth, number of reads supporting variant, ballence between forward and reverse reads

rule bcftools_filter_samtools:
    input:
        rules.anc_filter_samtools.output,
    output:
        f'{OUTPUT_DIR}/07_filtered/{{sample}}/{{sample}}_samtools_AB_AncFiltered.filt.vcf',
    conda:
        'envs/main.yml'
    shell:
        "bcftools filter -O v -o {output} -i 'MQ>30 & QUAL>75 & DP>40 & (DP4[2]+DP4[3])>4 & (DP4[0]+DP4[2])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.01 & (DP4[1]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.01' {input}"


rule bcftools_filter_samtools_two:
    input:
        rules.anc_filter_samtools.output,
    output:
        f'{OUTPUT_DIR}/07_filtered/{{sample}}/{{sample}}_samtools_filtered.vcf',
    conda:
        'envs/main.yml'
    shell:
        "bcftools filter -O v -o {output} -i 'MQ>30 & QUAL>75 & DP>10 & (DP4[2]+DP4[3])>4 & (DP4[2]+DP4[3])/DP>0.3 & (DP4[0]+DP4[2])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.01 & (DP4[1]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.01' {input}"


rule bcftools_filter_freebayes:
    input:
        rules.anc_filter_freebayes.output,
    output:
        f'{OUTPUT_DIR}/07_filtered/{{sample}}/{{sample}}_freebayes_BCBio_AncFiltered.filt.vcf',
    conda:
        'envs/main.yml'
    shell:
        "bcftools filter -O v -o {output} -i 'MQM>30 & MQMR>30 & QUAL>20 & INFO/DP>10 & (SAF+SAR)>4 & (SRF+SAF)/(INFO/DP)>0.01 & (SRR+SAR)/(INFO/DP)>0.01' {input}"


rule bcftools_filter_lofreq:
    input:
        rules.anc_filter_lofreq.output,
    output:
        f'{OUTPUT_DIR}/07_filtered/{{sample}}/{{sample}}_lofreq_tumor_relaxed_AncFiltered.filt.vcf',
    conda:
        'envs/main.yml'
    shell:
        "bcftools filter -O v -o {output} -i 'QUAL>20 & DP>20 & (DP4[2]+DP4[3])>4 & (DP4[0]+DP4[2])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.01 & (DP4[1]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.01' {input}"


#
# Fitler samtools by freebayes and lofreq
#
rule filter_samtools:
    input:
        samtools=rules.bcftools_filter_samtools.output,
        freebayes=rules.bcftools_filter_freebayes.output,
        lofreq=rules.bcftools_filter_lofreq.output,
    output:
        f'{OUTPUT_DIR}/07_filtered/{{sample}}/{{sample}}_samtools_AB_AncFiltered.filt.noOverlap.vcf',
    conda:
        'envs/main.yml'
    shell:
        "bedtools intersect -v -header -a {input.samtools} -b {input.freebayes} {input.lofreq} > {output}"


#
# Fitler freebayes by lofreq
#
rule filter_freebayes:
    input:
        freebayes=rules.bcftools_filter_freebayes.output,
        lofreq=rules.bcftools_filter_lofreq.output,
    output:
        f'{OUTPUT_DIR}/07_filtered/{{sample}}/{{sample}}_freebayes_BCBio_AncFiltered.filt.noOverlap.vcf',
    conda:
        'envs/main.yml'
    shell:
        "bedtools intersect -v -header -a {input.freebayes} -b {input.lofreq} > {output}"

#Filter by ancestor mutations
rule withoutanc_bcftools_filter_samtools_two:
	input:
	    before=rules.bcftools_filter_samtools_two.output,
	    anc=f'{OUTPUT_DIR}/07_filtered/anc/anc_samtools_filtered.vcf'
	output:
	    f'{OUTPUT_DIR}/07_filtered/{{sample}}/{{sample}}_samtools_filtered_noanc.vcf'
	conda:
	    'envs/main.yml'
	shell:
	    "bedtools intersect -v -header -a {input.before} -b {input.anc} > {output}"

rule withoutanc_bcftools_filter_lofreq:
	input:
	    before=rules.bcftools_filter_lofreq.output,
	    anc=f'{OUTPUT_DIR}/07_filtered/anc/anc_lofreq_tumor_relaxed_AncFiltered.filt.vcf'
	output:
	    f'{OUTPUT_DIR}/07_filtered/{{sample}}/{{sample}}_lofreq_tumor_relaxed_AncFiltered_noanc.filt.vcf'
	conda:
	    'envs/main.yml'
	shell:
	    "bedtools intersect -v -header -a {input.before} -b {input.anc} > {output}"

rule withoutanc_filter_samtools:
	input:
	    before=rules.filter_samtools.output,
	    anc=f'{OUTPUT_DIR}/07_filtered/anc/anc_samtools_AB_AncFiltered.filt.noOverlap.vcf'
	output:
	    f'{OUTPUT_DIR}/07_filtered/{{sample}}/{{sample}}_samtools_AB_AncFiltered_noanc.filt.noOverlap.vcf'
	conda:
	    'envs/main.yml'
	shell:
	    "bedtools intersect -v -header -a {input.before} -b {input.anc} > {output}"

rule withoutanc_filter_freebayes:
	input:
	    before=rules.filter_freebayes.output,
	    anc=f'{OUTPUT_DIR}/07_filtered/anc/anc_freebayes_BCBio_AncFiltered.filt.noOverlap.vcf'
	output:
	    f'{OUTPUT_DIR}/07_filtered/{{sample}}/{{sample}}_freebayes_BCBio_AncFiltered_noanc.filt.noOverlap.vcf'
	conda:
	    'envs/main.yml'
	shell:
	    "bedtools intersect -v -header -a {input.before} -b {input.anc} > {output}"


#Start annotating
rule annotate_samtools:
    input:
        rules.withoutanc_filter_samtools.output,
    output:
        f"{rules.filter_samtools.output}.annotated"
    conda:
        'envs/chris.yml'
    params:
        script=f'{workflow.basedir}/scripts/yeast_annotation_chris_edits_20170925.py'
    shell:
        f"python {{params.script}} -f {{input}} -s {config['annotate_sequences']} -n {config['annotate_noncoding']} -g {config['annotate_genome']}"


rule annotate_samtools_two:
    input:
        rules.withoutanc_bcftools_filter_samtools_two.output,
    output:
        f"{rules.bcftools_filter_samtools_two.output}.annotated"
    conda:
        'envs/chris.yml'
    params:
        script=f'{workflow.basedir}/scripts/yeast_annotation_chris_edits_20170925.py'
    shell:
        f"python {{params.script}} -f {{input}} -s {config['annotate_sequences']} -n {config['annotate_noncoding']} -g {config['annotate_genome']}"

rule annotate_freebayes:
    input:
        rules.withoutanc_filter_freebayes.output,
    output:
        f"{rules.filter_freebayes.output}.annotated"
    conda:
        'envs/chris.yml'
    params:
        script=f'{workflow.basedir}/scripts/yeast_annotation_chris_edits_20170925.py'
    shell:
        f"python {{params.script}} -f {{input}} -s {config['annotate_sequences']} -n {config['annotate_noncoding']} -g {config['annotate_genome']}"


rule annotate_lofreq:
    input:
        rules.withoutanc_bcftools_filter_lofreq.output,
    output:
        f"{rules.bcftools_filter_lofreq.output}.annotated"
    conda:
        'envs/chris.yml'
    params:
        script=f'{workflow.basedir}/scripts/yeast_annotation_chris_edits_20170925.py'
    shell:
        f"python {{params.script}} -f {{input}} -s {config['annotate_sequences']} -n {config['annotate_noncoding']} -g {config['annotate_genome']}"

rule lofreq_columns:
    input:
        rules.annotate_lofreq.output
    output:
        f"{rules.annotate_lofreq.output}".replace('.filt.vcf.annotated', '.filt.twoCol.vcf.annotated')
    conda:
        'envs/chris.yml'
    shell:
        r"""awk '$8 = $8 FS "NA NA"' """ + '{input} > {output}'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Run Pipeline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule finish:
    input:
        # QC steps
        rules.export_run_config.output,
        rules.list_samples.output,
        rules.run_fastqc_all.output,
        expand(rules.run_fastqc_persample.output, sample=SAMPLES),
        expand(rules.samtools_flagstat.output, sample=SAMPLES),
        # variant calling pipeline
        # KEEP realign.sort.bam.bai
        # KEEP realign.sort.bam
        # REMOVE EVERYTHING ELSE AT THE END, FIND A WAY TO DELETE EVERYTHING ELSE
        expand(rules.annotate_samtools.output, sample=SAMPLES),
        expand(rules.annotate_samtools_two.output, sample=SAMPLES),
        expand(rules.annotate_freebayes.output, sample=SAMPLES),
        expand(rules.lofreq_columns.output, sample=SAMPLES),
    output:
        f'{OUTPUT_DIR}/DONE.txt'
    shell:
        'touch {output}'



