#!/bin/bash

#SBATCH --array=1-68 #change according to number of sample to process
#SBATCH --job-name=leu_trim
#SBATCH --partition=cpu
#SBATCH --output=trim.%a.out
#SBATCH --error=trim.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --mail-user=sarah.schmid@unil.ch
#SBATCH --mail-type=ALL
#SBATCH --time=1-00:00:00

#we can change original name of file using the following command for example:
#ls -1 | cut -f 2 -d '-' | cut -f 1,3 -d _
#it will just keep the ID name and the R1 or R2 information


id=$(sed -n ${SLURM_ARRAY_TASK_ID}p samples.list) #e.g. for sample ID in samples.list: GB187

module load gcc/
module load trimmomatic/

trimmomatic PE /work/FAC/FBM/DBC/nsalamin/clownfish/sschmid/4_leucokranos/0_raw_reads/${id}_R1.fastq.gz /work/FAC/FBM/DBC/nsalamin/clownfish/sschmid/4_leucokranos/0_raw_reads/${id}_R2.fastq.gz ${id}_R1_paired.fq.gz ${id}_R1_unpaired.fq.gz ${id}_R2_paired.fq.gz ${id}_R2_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
