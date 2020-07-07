
#!/bin/bash
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=30G
#$ -l h_rt=24:00:00
#$ -l highp


########################### environment #####################################

. /u/home/c/christac/miniconda3/etc/profile.d/conda.sh
conda activate exome

mkdir -p gatk_output
mkdir -p gatk_output/samples

reference="../hg38/broad/Homo_sapiens_assembly38.fasta"
hapmap="/u/home/c/christac/project-nzaitlen/hg38/broad/hapmap_3.3.hg38.vcf.gz"
omni="/u/home/c/christac/project-nzaitlen/hg38/broad/1000G_omni2.5.hg38.vcf.gz"
thousand_genomes="/u/home/c/christac/project-nzaitlen/hg38/broad/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
dbsnp="/u/home/c/christac/project-nzaitlen/hg38/broad/Homo_sapiens_assembly38.dbsnp138.vcf"

indels="/u/home/c/christac/project-nzaitlen/hg38/broad/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

primary_targets="SeqCapEZ_Exome_v3.0_Design_Annotation_files/primary_targets.interval_list"

############################# combine samples  ##############################

# hard coded my samples here, idk a better way since each sample vcf needs to be a separate -V input life 
gatk GenomicsDBImport \
    -V gatk_output/41_variants.gvcf \
    -V gatk_output/42_variants.gvcf \
    -V gatk_output/43_variants.gvcf \
    -V gatk_output/44_variants.gvcf \
    -V gatk_output/45_variants.gvcf \
    -V gatk_output/46_variants.gvcf \
    -V gatk_output/47_variants.gvcf \
    -V gatk_output/48_variants.gvcf \
    -V gatk_output/49_variants.gvcf \
    -V gatk_output/50_variants.gvcf \
    -V gatk_output/51_variants.gvcf \
    --genomicsdb-workspace-path gatk_output/sample_database \
    --intervals intervals.list

gatk GenotypeGVCFs \
    -R $reference \
    -V gendb://gatk_output/sample_database \
    -O gatk_output/combined_samples.vcf 



############################ select variants ################################

gatk SelectVariants \
    -R $reference \
    -V gatk_output/combined_samples.vcf  \
    --select-type-to-include SNP \
    -O gatk_output/combined_samples_snps.vcf 


gatk SelectVariants \
    -R $reference \
    -V gatk_output/combined_samples.vcf  \
    --select-type-to-include INDEL \
    -O gatk_output/combined_samples_indels.vcf 


# ########################### recalibration ################################### 

 gatk VariantRecalibrator \
   -R $reference \
   -V gatk_output/combined_samples_snps.vcf \
   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
   -resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
   -resource:1000G,known=false,training=true,truth=false,prior=10.0 $thousand_genomes \
   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   -O gatk_output/combined_samples_snps.recal \
   --tranches-file gatk_output/combined_samples_snps.tranches \
   --rscript-file gatk_output/output.plots.R


gatk VariantRecalibrator \
    -R $reference \
    -V gatk_output/combined_samples_indels.vcf \
    -O gatk_output/combined_samples_indels.recal \
    --tranches-file gatk_output/combined_samples_indels.tranches \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 $indels \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode INDEL


# ######################## apply recalibration ##################################

 gatk ApplyVQSR \
   -R $reference \
   -V gatk_output/combined_samples_snps.vcf  \
   -O gatk_output/combined_samples_snps_recal.vcf \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file gatk_output/combined_samples_snps.tranches \
   --recal-file gatk_output/combined_samples_snps.recal \
   -mode SNP

gatk ApplyVQSR \
   -R $reference \
   -V gatk_output/combined_samples_indels.vcf  \
   -O gatk_output/combined_samples_indels_recal.vcf \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file gatk_output/combined_samples_indels.tranches \
   --recal-file gatk_output/combined_samples_indels.recal \
   -mode INDEL
 

# ############################ evaluate variants ##################################

gatk VariantEval \
    -R $reference \
    -L $primary_targets \
    -eval gatk_output/combined_samples_snps_recal.vcf \
    -D $dbsnp \
    -O gatk_output/combined_samples_snps.eval 

gatk VariantEval \
    -R $reference \
    -L $primary_targets \
    -eval gatk_output/combined_samples_indels_recal.vcf \
    -D $dbsnp \
    -O gatk_output/combined_samples_indels.eval 


# ######################## select variants  ##################################

gatk SelectVariants \
    -R $reference \
    -V gatk_output/combined_samples_snps_recal.vcf  \
    --select-type-to-include SNP \
    -O gatk_output/final_snps.vcf 

gatk SelectVariants \
    -R $reference \
    -V gatk_output/combined_samples_indels_recal.vcf  \
    --select-type-to-include INDEL \
    -O gatk_output/final_indels.vcf 







