#!/bin/bash

#SBATCH --job-name=mars_mafft_percula
#SBATCH --partition=cpu
#SBATCH --output=outscripts/%x.%a.out
#SBATCH --error=outscripts/%x.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=05:00:00

### This script is used to aligned the mitogenomes reconstructed with MITOBim
### The input file is the fasta file consisting in all the mitogenomes concatenated
#### Author: Sarah Schmid


##############################
# Arguments to pass
workindir=${1-/scratch/sschmi13/4_mitogenomes}	# working directory
input_fasta=${2-mitobim_outputs_percula/all_mtdna_percula.mitobim.fasta} # name of the input fasta file with all mitogenomes
outputdir=${3-aligned_mitogenomes}	# name of the output directory
output_prefix=${4-all_mtdna_percula}			# prefix name of the output files
thread=${5-4} #number of threads
mars=${6-./MARS/mars} #path to MARS binary executable

####################
# LOAD MODULES
####################

module load gcc/11.4.0 mafft/7.505

####################
# BEGINNNING OF SCRIPT
####################

# Create output directory
mkdir $workindir"/"$outputdir

#Replace x by n, because x is not recognised by MARS
sed -i 's/x/n/g' $input_fasta

# Running MARS
$mars -a DNA -i $input_fasta -o $workindir"/"$outputdir"/"$output_prefix".rotated.mars.fasta" -T 4

# Running mafft
mafft --thread 4 $workindir"/"$outputdir"/"$output_prefix".rotated.mars.fasta" > $workindir"/"$outputdir"/"$output_prefix".rotated.mars.mafft.fasta"