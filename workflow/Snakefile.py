
#
# yEvo Pipeline
#

from datetime import datetime

# create a new timestamped output directory for every pipeline run
OUTPUT_DIR = f"results/{datetime.now().strftime('%Y%m%d_%H%M%S')}"

# https://snakemake.readthedocs.io/en/v7.14.0/tutorial/basics.html#step-7-adding-a-target-rule 
rule all:
    input:
        f'{OUTPUT_DIR}/DONE.txt'




        


rule run_fastqc:
    input:
        config['fastq_dir']
    output:
        f'{OUTPUT_DIR}/01_fastqc/stdin_fastqc.html',
        f'{OUTPUT_DIR}/01_fastqc/stdin_fastqc.zip',
    conda:
        'envs/main.yml'
    shell:
        f'zcat {input}/*.fastq.gz | fastqc stdin --outdir={OUTPUT_DIR}/01_fastqc'
        

        
        
        
        

rule finish:
    input:
        rules.run_fastqc.output,
    output:
        f'{OUTPUT_DIR}/DONE.txt'
    conda:
        'envs/main.yml'
    shell:
        'touch {output}'
