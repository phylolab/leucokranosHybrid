#!/bin/bash

#SBATCH --array=1-66
#SBATCH --job-name=3_cla_processing
#SBATCH --partition=cpu
#SBATCH --output=processing.%a.out
#SBATCH --error=processing.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --mail-user=sarah.schmid@unil.ch
#SBATCH --mail-type=ALL
#SBATCH --time=1-00:00:00

id=$(sed -n ${SLURM_ARRAY_TASK_ID}p samples.list)

module load gcc
module load singularity
module load picard

export SINGULARITY_BINDPATH="/users,/work,/scratch"

python 3_processing.py /scratch/sschmi13/3_recomb_clarkii/3_processing/ ${id} /scratch/sschmi13/3_recomb_clarkii/2_mapping/mapping_sort_out/${id}.BWAMapped.Sorted.bam 4
#important to have the `/` at the end of the output path
