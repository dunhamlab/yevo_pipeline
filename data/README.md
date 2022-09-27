# Test Data for yEvo Pipeline

A complete set of pipeline inputs is included for testing the pipeline.
Provided files include required pipeline reference files (genome assembly, 
ancestral strain files, annotation files) and small example data files for 
three samples with 100 reads each:

```
data/
├── ancestor
│   ├── YMD4612_pink_S1_comb_R1R2.RG.MD.realign.sort.bam
│   ├── YMD4612_pink_S1_freebayes_BCBio.vcf
│   └── YMD4612_pink_S1_samtools_AB_AncFiltered.vcf
├── annotate
│   ├── orf_coding_all_R64-1-1_20110203.fasta
│   ├── S288C_reference_sequence_R64-1-1_20110203.fsa
│   └── saccharomyces_cerevisiae_R64-1-1_20110208.gff.filtered
├── fastq
│   ├── sample_01_R1_001.fastq.gz
│   ├── sample_01_R2_001.fastq.gz
│   ├── sample_02_R1_001.fastq.gz
│   ├── sample_02_R2_001.fastq.gz
│   ├── sample_03_R1_001.fastq.gz
│   └── sample_03_R2_001.fastq.gz
└── genome
    └── sacCer3.fasta
```

## Downloading test data

After following the installation instructions 
[here](https://github.com/dunhamlab/yevo_pipeline/blob/main/README.md#installation),
navigate to the repo's base directory and run:

```
$ ./scripts/download_data.sh
```
