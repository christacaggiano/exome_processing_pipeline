# Exome Processing Pipeline
A simple pipeline for processing exome sequencing data using BWA and GATK4. It takes raw fastq files and maps them to hg38 (or reference file of choice) and then performs basic QC. SNP and indels are called and processed using GATK4. While this pipeline is not completely automated (see usage instructions) it is lightweight and simple, and I developed it out of a need for a straightforward exome pipeline easily run on a linux distribution. 

I hope this pipeline is useful to you, but it certainly is not comprehensive- just what worked well for our application. Please email <christa@g.ucla.edu> with any suggestions or comments, or better yet, make a pull request! 

# Installation 

Requires `fastqc`, `multiqc`, `trimgalore`, `samtools`, `BWA`, `mosdepth`, `GATK4`, `R`, and `picard`. Required packages can be installed with the conda environment file `environment.yml` 

`conda env create -n exome -f environment.yml`

Note that the synatx between `GATK3` and `GATK4` differs, and amongst some versions of GATK4 (specifically in the `VariantRecalibrator` step). We use the version `gatk4-4.1.7.0-0` installed with conda.


# Use 

This pipeline was designed for a SGE cluster architecture, but can be easily adapted to any computing environment, as it is made of four relatively simple shell scripts. 

Assumes paired end fastq files in the input form: `$SAMPLE_NUMBER"_R1.fastq.gz"` and `$SAMPLE_NUMBER_"R2.fastq.gz"` The $SAMPLE_NUMBER parameter can be changed in the environment section of each script. 



## Preprocess 

Performs fastqc on raw fastq files, and then trims them.   

```bash 

# simple linux run (must be executable- use chmod a+x if not executable on your system)
./preproccess.sh


# SGE run
qsub preprocess.sh 

``` 


## Map 

Maps trimmed fastq files using BWA, allowing multimapping for downstream indel calling. This script will also calculate mapping statistics, along with on target read depth. This requires a bed file of target regions (ours is from [Roche](https://sequencing.roche.com/en/support-resources/discontinued-products/seqcap-ez-exome-v3-kit.html#:~:text=The%20SeqCap%C2%AE%20EZ%20Human,genes%20in%20the%20human%20genome)) and a reference genome fasta file. We use the newest [Broad hg38 version](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) downloaded on the command line using `gsutil`

We generate a histogram of target read depths using a custom R script `coverage_hist.R` and `tidyverse`. For histogram clarity, this script removes extreme read depths >1000x (extreme given an expected read depth approx 100x), but you might want to adapt given your needs. 

For speed, we run on a SGE HPC cluster with 8 threads, 10GB of memory per thread, and run all samples simulataneously. 


```bash 

# simple linux run (must be executable- use chmod a+x if not executable on your system)
./map.sh


# SGE run
qsub map.sh 

``` 

## Variant Identification 

This pipeline does variant identification in two steps. First, we call individual sample variants (SNPs and indels) using `GATK4 HaplotypeCaller`. These files are then combined into one vcf and post-processed together. 

To do this, we make use of several Broad reference files, again found [here].(https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1)

These files are:  

- hg38 reference: Homo_sapiens_assembly38.fasta
- hapmap 3: hapmap_3.3.hg38.vcf.gz
- 1000G omni: 1000G_omni2.5.hg38.vcf.gz
- 1000G SNPs: 1000G_phase1.snps.high_confidence.hg38.vcf.gz
- dbSNP: Homo_sapiens_assembly38.dbsnp138.vcf
- Mills and 1000G gold standard indels: Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
- hg38 known indels: Homo_sapiens_assembly38.known_indels.vcf.gz

The paths to these files can be specified at the beginning of each shell script. Please adapt the required files according to your needs. 

Here, the target file discussed in the `map.sh` section needs to be in GATK's interval list format. To do this, we used `Picard BedToIntervalList`. This required both a bed file of our primary targets and a hg38 dictionary file, found in the Broad Resource Bundle (see above). 

The requisite command is roughly: 
```bash
picard BedToIntervalList \ 
 I=SeqCap_EZ_Exome_v3_hg38_primary_targets.bed \ 
 O=SeqCap_EZ_Exome_v3_hg38_primary_targets.interval_list \ 
 SD=../../hg38/broad/Homo_sapiens_assembly38.dic 
```


### Step 1: identification 

The GATK4 applications are quite memory intensive. We use 30G of memory (single threaded) for each sample, run conconcurrently. This may need to be adapted based on your starting bam file size. 

```bash 

# simple linux run (commands are very memory intensive, so pipeline may fail depending on available memory)
./call_variants.sh


# SGE run
qsub call_variants.sh 

``` 

### Step 2: final variants 

This step merges all sample gvcfs into one file, as recommended by the GATK developers. It then uses `GATK4 VariantRecalibrator` which is a machine learning algorithm that learns systematic biases in the data and corrects them. The more samples, the better. 

Currently, the command to merge all gvcfs together is hardcoded with the sample name, so please change the`GenomicsDBImport` command with the correct file names (or pull request if you have a better way of doing this)

Again, we use 30GB of memory, but since we are combining all samples, we only run one command. 

```bash 

# simple linux run (commands are very memory intensive, so pipeline may fail depending on available memory)
./finalize_variants.sh


# SGE run
qsub finalize_variants.sh 

``` 


