#!/bin/bash

#SBATCH --array=1-24
#SBATCH --job-name=plot_individual_elai
#SBATCH --partition=cpu
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G
#SBATCH --time=02:00:00

# Set the path for R to find and load the libraries

module load gcc/11.4.0 r/4.3.2

#export R_LIBS_USER=/users/sschmi13/R/x86_64-pc-linux-gnu-library/4.3 # The needed packages should be already installed in this directory

id=$(sed -n ${SLURM_ARRAY_TASK_ID}p chr.list)

Rscript 6_${id}_elai_individual_plot.R