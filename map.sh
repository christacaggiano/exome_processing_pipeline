
#!/bin/bash
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=10G
#$ -l h_rt=24:00:00
#$ -l highp
#$ -t 1-11
#$ -pe shared 8
######
######### #$ -pe shared 4

########################### environment #####################################

. /u/home/c/christac/miniconda3/etc/profile.d/conda.sh
conda activate exome

ID_NUM=$((SGE_TASK_ID+40))
mkdir -p aligned
mkdir -p aligned/coverage
mkdir -p trash


reference="../hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"  # reference fasta file (I got mine from https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1)
primary_targets="SeqCapEZ_Exome_v3.0_Design_Annotation_files/primary_targets.bed"  # exome target file in bed format

########################### mapping #####################################
## index genome - only need to do this step once 
bwa index $reference -a bwtsw

## map to hg38 genome with alt contigs using 8 threads
bwa mem -t 8 $reference  "trimmed/"$ID_NUM"_R1_val_1.fq.gz" "trimmed/"$ID_NUM"_R2_val_2.fq.gz" > "aligned/"$ID_NUM".sam"

## convert sam to bam and sort 
samtools view -S -b "aligned/"$ID_NUM".sam" > "aligned/"$ID_NUM".bam"
samtools collate -@ 4  "aligned/"$ID_NUM".bam" "aligned/"$ID_NUM".collated"
samtools fixmate -m -@ 4 "aligned/"$ID_NUM".collated.bam" "aligned/"$ID_NUM".fixmate.bam" 
samtools index "aligned/"$ID_NUM".nodup.bam"

########################### sort and mark duplicates #####################################

samtools sort -@ 4 "aligned/"$ID_NUM".fixmate.bam" -o "aligned/"$ID_NUM".fixmate.sorted.bam"

## picard markdups
picard MarkDuplicates I="aligned/"$ID_NUM".fixmate.sorted.bam"  O="aligned/"$ID_NUM".marked_duplicates.bam" M="aligned/"$ID_NUM".marked_dup_metrics.txt"

# samtools sort -n -@ 4 "aligned/"$ID_NUM"marked_duplicates.bam" -o "aligned/"$ID_NUM".sorted.bam"
samtools index "aligned/"$ID_NUM".marked_duplicates.bam"

picard AddOrReplaceReadGroups \
       I="aligned/"$ID_NUM".marked_duplicates.bam" \
       O="aligned/"$ID_NUM".readgroup.bam" \
       RGID=$ID_NUM \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=$ID_NUM
 
 picard ValidateSamFile \
      I="aligned/"$ID_NUM".readgroup.bam" \
      MODE=SUMMARY

samtools index "aligned/"$ID_NUM".readgroup.bam"

########################### cleanup  #####################################
# remove intermediate files (For safety, I moved mine to a trash folder, which can later be deleted)
mv "aligned/"$ID_NUM".sam" trash
mv "aligned/"$ID_NUM".collated" trash
mv "aligned/"$ID_NUM".fixmate.bam" trash
mv "aligned/"$ID_NUM".nodup.bam" trash 

########################### QC  #####################################

## mapping percentage
echo ""
echo "fastq reads"
echo $ID_NUM"_R1.fq.gz"
echo $(zcat $ID_NUM"_R1.fq.gz"|wc -l)/4|bc

echo ""
echo "trimmed reads"
echo "trimmed/"$ID_NUM"_R1_val_1.fq.gz"
echo $(zcat "trimmed/"$ID_NUM"_R1_val_1.fq.gz"|wc -l)/4|bc

echo ""
echo "aligned reads"
samtools flagstat "aligned/"$ID_NUM".readgroup.bam"

# unmapped reads 
samtools view -b -f 4 aligned/$ID_NUM".sam" > aligned/$ID_NUM".unmapped.bam"
samtools flagstat aligned/$ID_NUM".unmapped.bam"


################################## some more QC #####################################

# picard summary stats
picard \
        CollectAlignmentSummaryMetrics \
        R=$reference  \
        I="aligned/"$ID_NUM".readgroup.bam" \
        O="aligned/coverage/"$ID_NUM".alignment_metrics.txt"


# insert size statistics 
picard \
CollectInsertSizeMetrics \
        INPUT=s"aligned/"$ID_NUM".readgroup.bam" \
        OUTPUT="aligned/coverage/"$ID_NUM".insert_metrics.txt"\
        HISTOGRAM_FILE="aligned/coverage/"$ID_NUM".insert_size_hist.pdf"


## read depth at targets
mosdepth -t 4 --by $primary_targets "aligned/coverage/"$ID_NUM"_mosdepth" "aligned/"$ID_NUM".readgroup.bam"

gunzip "aligned/coverage/"$ID_NUM"_mosdepth.regions.bed.gz"

# plots target read depth file
Rscript --vanilla coverage_hist.R $ID_NUM "aligned/coverage/"$ID_NUM"_mosdepth.regions.bed"


