
#!/bin/bash
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=15G
#$ -l h_rt=6:00:00
#$ -l highp
#$ -t 


. /u/home/c/christac/miniconda3/etc/profile.d/conda.sh
conda activate exome

ID_NUM=$((SGE_TASK_ID+40))
raw_data_folder="raw_data"


########################### QC ##############################

# fastqc on individual fastq files 
fastqc $raw_data_folder"/"$ID_NUM*".fq.gz" --outdir='fastqc/' 
multiqc 'fastqc/' 

# trimgalore on combined fastq files 
trim_galore --paired "fastqc/"$ID_NUM"_R1.fq.gz" $ID_NUM"_R2.fq.gz" 

# fastqc on trimmed files 
fastqc "fastqc/"$ID_NUM*"_val_"*".fq.gz" --outdir='fastqc_trimmed/' 
multiqc 'fastqc_trimmed/'