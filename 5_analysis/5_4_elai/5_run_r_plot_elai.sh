#!/usr/bin/bash -l

#SBATCH --job-name=plot_elai
#SBATCH --partition=cpu
#SBATCH --output=outscripts/%x.out
#SBATCH --error=outscripts/%x.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=03:00:00

# Set the path for R to find and load the libraries
# export R_LIBS_USER=/users/sschmi13/R/x86_64-pc-linux-gnu-library/4.3 # The needed packages should be already installed in this directory

module load gcc/11.4.0 r/4.3.2

Rscript Scripts/ELAI_allchr_ggbio.R

