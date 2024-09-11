#!/usr/bin/bash -l

#SBATCH --job-name=BWA
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --output=outscripts/BWAIndex.out
#SBATCH --error=outscripts/BWAIndex.err
#SBATCH --time=1:00:00

# Use this script only if the Reference indexing (for BWA) has not been performed yet

##############################
# Argument to pass: 
pathtoReference=${1-0_AclarkiiReference}
RefGenome=${2-AclarkiiGenome.Chr.fna} # GCA_027123335.1_OIST_Acla_v1_genomic.fna.gz was downloaded from 
							# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/123/335/GCA_027123335.1_OIST_Acla_v1
							# and was renamed as AclarkiiGenome.fna.gz. AclarkiiGenome.Chr.fna contains only
							# assembled chromosomes with modified names (e.g. Chr1, Chr2...etc)
RefOutputName=${3-AclarkiiGenome.Chr}
##############################

# Check settings
echo "Using following input:
Path to reference genome:		"$pathtoReference"
Reference Genome Name:			"$RefGenome"
Reference output name for indexing: 	"$RefOutputName

echo "If settings are wrong, please look at instructions

Indexing..."

####################
# LOAD MODULES
####################
module load bwa/0.7.17


####################
# BEGINNNING OF SCRIPT
####################
bwa index $pathtoReference"/"$RefGenome -p $pathtoReference"/"$RefOutputName

####################
# END
####################


echo "Finish!"
