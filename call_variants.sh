
#!/bin/bash
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=30G
#$ -l h_rt=24:00:00
#$ -l highp
#$ -t 1-11


########################### environment #####################################

. /u/home/c/christac/miniconda3/etc/profile.d/conda.sh
conda activate exome

ID_NUM=$((SGE_TASK_ID+40))
mkdir -p gatk_output
mkdir -p gatk_output/bams

reference="../hg38/broad/Homo_sapiens_assembly38.fasta"
hapmap="/u/home/c/christac/project-nzaitlen/hg38/broad/hapmap_3.3.hg38.vcf.gz"
omni="/u/home/c/christac/project-nzaitlen/hg38/broad/1000G_omni2.5.hg38.vcf.gz"
thousand_genomes="/u/home/c/christac/project-nzaitlen/hg38/broad/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
dbsnp="/u/home/c/christac/project-nzaitlen/hg38/broad/Homo_sapiens_assembly38.dbsnp138.vcf"

mills="/u/home/c/christac/project-nzaitlen/hg38/broad/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
indels="/u/home/c/christac/project-nzaitlen/hg38/broad/Homo_sapiens_assembly38.known_indels.vcf.gz"

primary_targets="SeqCapEZ_Exome_v3.0_Design_Annotation_files/primary_targets.interval_list"


# ########################### base recalibration  #####################################

# recalibrates for indel identification 
gatk BaseRecalibrator \
        -R $reference \
        -I "aligned/"$ID_NUM".readgroup.bam" \
        -O "gatk_output/bams/"$ID_NUM"_recal.table" \
        --known-sites $indels \
        --known-sites $mills \
        --known-sites $dbsnp 

gatk ApplyBQSR \
        -R $reference \
        -I "aligned/"$ID_NUM".readgroup.bam" \
        -bqsr "gatk_output/bams/"$ID_NUM"_recal.table" \
        -O "gatk_output/bams/"$ID_NUM"_recal_reads.bam" \


gatk BaseRecalibrator \
        -R $reference \
        -I "gatk_output/bams/"$ID_NUM"_recal_reads.bam"  \
        -O "gatk_output/bams/"$ID_NUM"_post_recal.table" \
        --known-sites $indels \
        --known-sites $mills \
        --known-sites $dbsnp 

# plots that help to analyze success of base recalibration 
gatk AnalyzeCovariates \
    -before "gatk_output/bams/"$ID_NUM"_recal.table" \
    -after "gatk_output/bams/"$ID_NUM"_post_recal.table" \
    -plots "gatk_output/bams/"$ID_NUM"_recalibration_plots.pdf"


############################# variant calling #########################################

# calls sample snps and indels in gvcf form 

gatk HaplotypeCaller \
        -R $reference \
        -D $dbsnp \
        -I "gatk_output/bams/"$ID_NUM"_recal_reads.bam" \
        -ERC GVCF \
        -ip 100 \
        -L $primary_targets \
        -O "gatk_output/"$ID_NUM"_variants.gvcf"






