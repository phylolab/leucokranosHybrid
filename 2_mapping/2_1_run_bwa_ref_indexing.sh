#!/bin/bash -l

#SBATCH --job-name=BWA
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --output=outscripts/bwa_index.out
#SBATCH --error=outscripts/bwa_index.err
#SBATCH --time=1:00:00

# Use this script only if the reference indexing (for BWA) has not been performed yet
# Author: Wan-Ting Huang

##############################
# Arguments to pass 
path_to_ref=${1-0_aclarkii_ref}
ref_output_name=${3-aclarkii_genome.chr}
ref_genome=${2-aclarkii_genome.chr.fna} 
# GCA_027123335.1_OIST_Acla_v1_genomic.fna.gz was downloaded from 
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/123/335/GCA_027123335.1_OIST_Acla_v1
# and was renamed as AclarkiiGenome.fna.gz. aclarkii_genome.chr.fna contains only
# assembled chromosomes with modified names (e.g. Chr1, Chr2...etc)
##############################

# Check settings
echo "Using following input:
Path to reference genome:		"$path_to_ref"
Reference Genome Name:			"$ref_genome"
Reference output name for indexing: 	"$ref_output_name

echo "If settings are wrong, please look at instructions

Indexing..."

####################
# LOAD MODULES
####################
module load bwa/0.7.17


####################
# BEGINNNING OF SCRIPT
####################
bwa index $path_to_ref"/"$ref_genome -p $path_to_ref"/"$ref_output_name

####################
# END
####################


echo "Finish! Please run 2_2_run_bwa_mapping.sh to map reads against reference genome"
