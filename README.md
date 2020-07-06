# Exome Processing Pipeline
A simple pipeline for processing exome sequencing data using BWA and GATK4. Takes raw fastq files and maps them to hg38 (or reference file of choice) and then performs basic QC. SNP and indels are called and processed using GATK4. 

# Installation 

Requires `fastqc`, `multiqc`, `trimgalore`, `samtools`, `BWA`, `GATK4`, `R`, and `picard`. Required packages can be installed with the conda environment file `environment.yml` 

`conda env create -n exome -f environment.yml`


# Use 

This pipeline was designed for a SGE cluster architecture, but can be easily adapted to any computing environment, as it is made of four relatively simple shell scripts. 

Assumes paired end fastq files in the input form: `$SAMPLE_NUMBER"_R1.fastq.gz"` and `$SAMPLE_NUMBER_"R2.fastq.gz"` The $SAMPLE_NUMBER parameter can be changed in the environment section of each script. 



## Preprocess 

Performs fastqc on raw fastq files, and then trims them.   
```bash 

# simple linux run (must be executable- use chmod a+x if not executable on your system)
./preproccess


# SGE run
qsub preprocess.sh 
``` 


