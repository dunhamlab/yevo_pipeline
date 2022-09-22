
# gatk-3.7 conda installation notes

### The problem

After creating the conda env with the `gatk=3.7` package from bioconda, GATK still cannot be used:

```
$ gatk
GATK jar file not found. Have you run "gatk-register"?
```

Running the provided command produces:

```
$ gatk-register
  It looks like GATK has not yet been installed.

  Usage: gatk-register /path/to/GenomeAnalysisTK[-3.7.tar.bz2|.jar]

  Due to license restrictions, this recipe cannot distribute 
  and install GATK directly. To fully install GATK, you must 
  download a licensed copy of GATK from the Broad Institute: 
    https://www.broadinstitute.org/gatk/download/ 
  and run (after installing this package):
    gatk-register /path/to/GenomeAnalysisTK[-3.7.tar.bz2|.jar], 
   This will copy GATK into your conda environment.
```


### The solution

1. The source file `GenomeAnalysisTK-3.7-0-gcfedb67.tar.bz2` was manually downloaded (on 2022-09-21) from https://console.cloud.google.com/storage/browser/_details/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.7-0-gcfedb67.tar.bz2

2. The source file is included in this repo 

3. Snakemake runs the register command (in the `register_gatk` rule) on each pipeline run (redundantly, so that there's no one-time environment setup needed upon fresh repo clones etc.)


**NOTE:** This is a temp fix; alternative solutions where the running pipeline doesn't attempt to modify its existing conda env should be explored (potential issues with concurrent pipeline runs? potential solutions with other complete versions of gatk or custom conda build?)