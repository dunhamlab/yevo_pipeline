#
# align_compare.py
#
# authors: Lucas Kampman, ChatGPT=4.0
#
# date: 2023/09/21
#
# Usage:
# python align_compare.py --ref ref.fasta \
#                           --ancestral_1 ancestral_1.fastq \
#                           --ancestral_2 ancestral_2.fastq \
#                           --csv evolved_samples.csv \
#                           [--indexed] 
# The --indexed flag can be used if the indexed genome already 
# exists in the working directory.
#
#
# example evolved_samples.csv:
# (NOTE: this file is expected to have a header line.)
#
# sample_name,paired-end-read1,paired-end-read2 
# evolved1,data/evolved1_R1.fastq,data/evolved1_R2.fastq
# evolved2,data/evolved2_R1.fastq,data/evolved2_R2.fastq
# etc.
#
import subprocess
import csv
def run_cmd(cmd):
    """Execute the command."""
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    if process.returncode != 0:
        print(f"Error executing: {cmd}")
        print(err.decode())
        exit(1)
    return out
def index_ref_genome(ref):
    # Index the reference genome
    run_cmd(f"bwa index {ref}")
    run_cmd(f"samtools faidx {ref}")
def align_ancestral(ref, ancestral_1, ancestral_2):
    # Align paired-end reads for ancestral sample
    run_cmd(f"bwa mem {ref} {ancestral_1} {ancestral_2} > ancestral.sam")
    # Convert SAM to BAM and sort
    run_cmd("samtools view -Sb ancestral.sam | samtools sort -o ancestral_sorted.bam")
    # Index the BAM file
    run_cmd("samtools index ancestral_sorted.bam")
def align_evolved(ref, sample_name, evolved_1, evolved_2):
    # Align paired-end reads for evolved sample
    run_cmd(f"bwa mem {ref} {evolved_1} {evolved_2} > {sample_name}.sam")
    # Convert SAM to BAM and sort
    run_cmd(f"samtools view -Sb {sample_name}.sam | samtools sort -o {sample_name}_sorted.bam")
    # Index BAM files
    run_cmd(f"samtools index {sample_name}_sorted.bam")
def call_variants_bcftools(ref, sample_name):
    # Variant calling using BCFtools
    run_cmd(f"bcftools mpileup -O b -o {sample_name}_bcftools_raw.bcf -f {ref} {sample_name}_sorted.bam")
    run_cmd(f"bcftools call -vmO z -o {sample_name}_bcftools_filtered.vcf.gz {sample_name}_bcftools_raw.bcf")
def call_variants_gatk(ref, sample_name):
    # Assuming GATK is installed and the "gatk" command is in your PATH
    # Variant calling using GATK's HaplotypeCaller
    run_cmd(f"gatk HaplotypeCaller -R {ref} -I {sample_name}_sorted.bam -O {sample_name}_GATK_filtered.vcf.gz")
def call_variants_freebayes(ref, sample_name):
    # Variant calling using FreeBayes
    run_cmd(f"freebayes -f {ref} {sample_name}_sorted.bam > {sample_name}_FreeBayes_filtered.vcf")
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Align paired-end reads to a reference genome and call variants.")
    parser.add_argument("--ref", required=True, help="Path to the reference genome fasta file.")
    parser.add_argument("--ancestral_1", required=True, help="Path to the ancestral paired-end reads (file 1).")
    parser.add_argument("--ancestral_2", required=True, help="Path to the ancestral paired-end reads (file 2).")
    parser.add_argument("--csv", required=True, help="Path to the CSV file with evolved sample details.")
    parser.add_argument("--indexed", action='store_true', help="Specify this flag if the reference genome is already indexed.")
    args = parser.parse_args()
    # If genome is not indexed, index it
    if not args.indexed:
        index_ref_genome(args.ref)
    # Align ancestral
    align_ancestral(args.ref, args.ancestral_1, args.ancestral_2)
    # Loop over evolved samples, align them, and call variants using multiple methods
    with open(args.csv, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        header = next(csv_reader)  # Skip header
        for row in csv_reader:
            sample_name, evolved_1, evolved_2 = row
            align_evolved(args.ref, sample_name, evolved_1, evolved_2)
            call_variants_bcftools(args.ref, sample_name)
            #call_variants_gatk(args.ref, sample_name)
            #call_variants_freebayes(args.ref, sample_name)
