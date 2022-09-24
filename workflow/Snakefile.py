
#
# yEvo Pipeline
#

import os
from datetime import datetime

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Define Constants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# discover input files
SAMPLES = list(set(glob_wildcards(f"{config['fastq_dir']}/{{sample}}_R1_001.fastq.gz").sample))

# create a new timestamped output directory for every pipeline run
OUTPUT_DIR = f"results/{datetime.now().strftime('%Y%m%d_%H%M%S')}"
# OUTPUT_DIR = f"results/20220923_024004"

# Project name and date for bam header
SEQID='yevo_pipeline_align'


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Begin Pipeline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# https://snakemake.readthedocs.io/en/v7.14.0/tutorial/basics.html#step-7-adding-a-target-rule 
rule all:
    input:
        f'{OUTPUT_DIR}/DONE.txt'


# ~~~~~~~~~~~~~~~~~~~~~~~~~~ Set Up Reference Files ~~~~~~~~~~~~~~~~~~~~~~~~~ #

#
# make a list of discovered samples
#
rule list_samples:
    output:
        f"{OUTPUT_DIR}/01_ref_files/00_sample_list.txt"
    shell:
        "echo -e '{}' > {{output}}".format('\n'.join(SAMPLES))


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


rule index_fasta:
    input:
        rules.copy_fasta.output
    output:
        f"{rules.copy_fasta.output}.fai"
    conda:
        'envs/main.yml'
    shell:
        "samtools faidx {input}"


rule create_ref_dict:
    input:
        rules.copy_fasta.output
    output:
        f"{rules.copy_fasta.output}".rstrip('fasta') + 'dict'
    conda:
        'envs/main.yml'
    shell:
        "picard CreateSequenceDictionary -R {input}"

#
# copy the supplied reference genome fasta to the pipeline output directory for reference
#
rule copy_ancestor_bam:
    input:
        config['ancestor_bam']
    output:
        f"{OUTPUT_DIR}/01_ref_files/{os.path.basename(config['ancestor_bam'])}"
    shell:
        "cp {input} {output}"


rule index_ancestor_bam:
    input:
        rules.copy_ancestor_bam.output
    output:
        f"{rules.copy_ancestor_bam.output}.bai"
    conda:
        'envs/main.yml'
    shell:
        "samtools index {input}"


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
        r"""bwa mem -R '@RG\tID:""" + SEQID + r"""\tSM:""" + '{wildcards.sample}' + r"""\tLB:1'""" + ' {rules.copy_fasta.output} {input.r1} {input.r2} > {output}'


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
        bam=rules.samtools_sort_one.output,
        idx=rules.samtools_index_one.output
    output:
        f"{OUTPUT_DIR}/03_init_alignment/{{sample}}_R1R2_sort_flagstat.txt"
    conda:
        'envs/main.yml'
    shell:
        "samtools flagstat {input.bam} > {output}"


rule picard_mark_dupes:
    input:
        bam=rules.samtools_sort_one.output,
        idx=rules.samtools_index_one.output,
        dct=rules.create_ref_dict.output
    output:
        bam=f"{OUTPUT_DIR}/04_picard/{{sample}}_comb_R1R2.MD.bam",
        metrics=f"{OUTPUT_DIR}/04_picard/{{sample}}_comb_R1R2.sort_dup_metrix"
    conda:
        'envs/main.yml'
    shell:
        "picard MarkDuplicates --INPUT {input.bam} --OUTPUT {output.bam} --METRICS_FILE {output.metrics} --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY LENIENT"


rule picard_read_groups:
    input:
        rules.picard_mark_dupes.output.bam
    output:
        f"{OUTPUT_DIR}/04_picard/{{sample}}_comb_R1R2.RG.MD.bam",
    conda:
        'envs/main.yml'
    shell:
        f"picard AddOrReplaceReadGroups --INPUT {{input}} --OUTPUT {{output}} --RGID {SEQID} --RGLB 1 --RGPU 1 --RGPL illumina --RGSM {{wildcards.sample}} --VALIDATION_STRINGENCY LENIENT"



rule samtools_sort_two:
    input:
        rules.picard_read_groups.output
    output:
        f"{OUTPUT_DIR}/04_picard/{{sample}}_comb_R1R2.RG.MD.sort.bam",
    conda:
        'envs/main.yml'
    shell:
        "samtools sort {input} -o {output}"


rule samtools_index_two:
    input:
        rules.samtools_sort_two.output
    output:
        f"{OUTPUT_DIR}/04_picard/{{sample}}_comb_R1R2.RG.MD.sort.bam.bai",
    conda:
        'envs/main.yml'
    shell:
        "samtools index {input}"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GATK Re-alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#
#  configure gatk-3.7 inside conda environment
#
rule gatk_register:
    input:
        f'workflow/envs/src/GenomeAnalysisTK-3.7-0-gcfedb67.tar.bz2'
    output:
        f'{OUTPUT_DIR}/05_gatk/gatk_3.7_registered.txt'
    conda:
        'envs/main.yml'
    shell:
        "gatk-register {input} > {output}"


rule gatk_realign_targets:
    input:
        fa=rules.copy_fasta.output,
        bam=rules.samtools_sort_two.output,
        idx=rules.samtools_index_two.output,
        gatk=rules.gatk_register.output,
        faidx=rules.index_fasta.output
    output:
        f'{OUTPUT_DIR}/05_gatk/{{sample}}_comb_R1R2.bam.intervals'
    conda:
        'envs/main.yml'
    shell:
        "GenomeAnalysisTK -T RealignerTargetCreator -R {input.fa} -I {input.bam} -o {output}"


rule gatk_realign_indels:
    input:
        fa=rules.copy_fasta.output,
        intervals=rules.gatk_realign_targets.output,
        bam=rules.samtools_sort_two.output,
        idx=rules.samtools_index_two.output        
    output:
        bam=f'{OUTPUT_DIR}/05_gatk/{{sample}}_comb_R1R2.RG.MD.realign.bam',
        bai=f'{OUTPUT_DIR}/05_gatk/{{sample}}_comb_R1R2.RG.MD.realign.bai'
    conda:
        'envs/main.yml'
    shell:
        "GenomeAnalysisTK -T IndelRealigner -R {input.fa} -I {input.bam} -targetIntervals {input.intervals} -o {output.bam}"


rule samtools_sort_three:
    input:
        rules.gatk_realign_indels.output.bam
    output:
        f'{OUTPUT_DIR}/05_gatk/{{sample}}_comb_R1R2.RG.MD.realign.sort.bam',
    conda:
        'envs/main.yml'
    shell:
        "samtools sort {input} -o {output}"


rule samtools_index_three:
    input:
        rules.samtools_sort_three.output
    output:
        f'{OUTPUT_DIR}/05_gatk/{{sample}}_comb_R1R2.RG.MD.realign.sort.bam.bai',
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
        f'{OUTPUT_DIR}/06_variant_calling/{{sample}}_samtools_AB.vcf',
    conda:
        'envs/main.yml'
    shell:
        "bcftools mpileup --ignore-RG -Ou -ABf {rules.copy_fasta.output} {input.bam} | bcftools call -vmO v -o {output}"


rule freebayes:
    input:
        bam=rules.samtools_sort_three.output,
        idx=rules.samtools_index_three.output
    output:
        f'{OUTPUT_DIR}/06_variant_calling/{{sample}}_freebayes_BCBio.vcf',
    conda:
        'envs/main.yml'
    shell:
        "freebayes -f {rules.copy_fasta.output} --pooled-discrete --pooled-continuous --report-genotype-likelihood-max --allele-balance-priors-off --min-alternate-fraction 0.1 {input.bam} > {output}"


rule lofreq:
    input:
        bam=rules.samtools_sort_three.output,
        idx=rules.samtools_index_three.output,
        ancidx=rules.index_ancestor_bam.output,
    output:
        normal=f'{OUTPUT_DIR}/06_variant_calling/{{sample}}_lofreq_normal_relaxed.vcf.gz',
        tumor=f'{OUTPUT_DIR}/06_variant_calling/{{sample}}_lofreq_tumor_relaxed.vcf.gz',
        somatic=f'{OUTPUT_DIR}/06_variant_calling/{{sample}}_lofreq_somatic_final.snvs.vcf.gz',
    conda:
        'envs/main.yml'
    shell:
        f"lofreq somatic -n {{rules.copy_ancestor_bam.output}} -t {{input.bam}} -f {{rules.copy_fasta.output}} -o {OUTPUT_DIR}/06_variant_calling/{{wildcards.sample}}_lofreq_"


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
        rules.bcftools_pileup.output,
    output:
        f'{OUTPUT_DIR}/07_filtered/{{sample}}_samtools_AB_AncFiltered.vcf',
    conda:
        'envs/main.yml'
    shell:
        f"bedtools intersect -v -header -a {{input}} -b {config['ancestor_samtools_vcf']} > {{output}}"


#
# filter freebayes results by ancestor
#
rule anc_filter_freebayes:
    input:
        rules.freebayes.output,
    output:
        f'{OUTPUT_DIR}/07_filtered/{{sample}}_freebayes_BCBio_AncFiltered.vcf',
    conda:
        'envs/main.yml'
    shell:
        f"bedtools intersect -v -header -a {{input}} -b {config['ancestor_freebayes_vcf']} > {{output}}"


#
# filter lofreq results by ancestor
#
rule anc_filter_lofreq:
    input:
        normal=rules.unzip_lofreq.output.normal,
        tumor=rules.unzip_lofreq.output.tumor,
    output:
        f'{OUTPUT_DIR}/07_filtered/{{sample}}_lofreq_tumor_relaxed_AncFiltered.vcf',
    conda:
        'envs/main.yml'
    shell:
        f"bedtools intersect -v -header -a {{input.tumor}} -b {{input.normal}} > {{output}}"










rule finish:
    input:
        rules.list_samples.output,
        rules.run_fastqc_all.output,
        expand(rules.run_fastqc_persample.output, sample=SAMPLES),
        
        expand(rules.samtools_flagstat.output, sample=SAMPLES),
        
        expand(rules.anc_filter_samtools.output, sample=SAMPLES),
        expand(rules.anc_filter_freebayes.output, sample=SAMPLES),
        expand(rules.anc_filter_lofreq.output, sample=SAMPLES),
        

    output:
        f'{OUTPUT_DIR}/DONE.txt'
    shell:
        'touch {output}'
