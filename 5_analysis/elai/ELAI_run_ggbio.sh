#!/usr/bin/bash -l

#SBATCH --job-name=allchr_ggbio
#SBATCH --partition=cpu
#SBATCH --output=outscripts/allchr_ggbio.out
#SBATCH --error=outscripts/allchr_ggbio.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=03:00:00

# Set the path for R to find and load the libraries
export R_LIBS_USER=/scratch/whuang4/RLIB_elai # The needed packages should be already installed in this directory

module load r-light/4.4.1

Rscript Scripts/ELAI_allchr_ggbio.R