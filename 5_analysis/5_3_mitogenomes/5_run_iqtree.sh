#!/bin/bash

#SBATCH --job-name=iqtree
#SBATCH --partition=cpu
#SBATCH --output=outscripts/%x.%A.out
#SBATCH --error=outscripts/%x.%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=1-00:00:00

### This script is used to reconstruct phylogenies based on the mitogenomes reconstructed with MITObim and aligned with mafft
### The input file is the fasta file consisting in all the mitogenomes concatenated and aligned
#### Author: Sarah Schmid


##############################
# Arguments to pass
workindir=${1-/scratch/sschmi13/4_mitogenomes}	# working directory
fasta=${2-aligned_mitogenomes/all_mtdna_percula.rotated.mars.mafft.fasta} # alignment in fasta
boot=${3-1000} # bootstraps
output=${4-iqtree_output/all_mtdna_percula} #output dir


module load gcc/11.4.0 iqtree2/2.2.2

iqtree2 -s $fasta -alrt 1000 -b $boot --prefix $output

