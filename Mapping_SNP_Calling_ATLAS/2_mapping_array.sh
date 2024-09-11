#!/bin/bash

#SBATCH --array=1-30
#SBATCH --job-name=2_mapping
#SBATCH --partition=cpu
#SBATCH --output=mapping.%a.out
#SBATCH --error=mapping.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --mail-user=sarah.schmid@unil.ch
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00

id=$(sed -n ${SLURM_ARRAY_TASK_ID}p samples.list)

module load gcc/
module load bwa/

mkdir mapping_out

python 2_mapping.py mapping_out ${id} /work/FAC/FBM/DBC/nsalamin/clownfish/sschmid/reference_genome/Amphiprion_percula.Nemo_v1.dna.toplevel.fa /scratch/sschmi13/3_recomb_clarkii/1_trimming/${id}_R1_paired.fq.gz /scratch/sschmi13/3_recomb_clarkii/1_trimming/${id}_R2_paired.fq.gz 4
