#!/usr/bin/bash -l

#SBATCH --job-name=index_ref_gatk
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --output=outscripts/index_ref_gatk.out
#SBATCH --error=outscripts/index_ref_gatk.err
#SBATCH --time=01:00:00

#### The script is used to generate index files for performing GATK HaplotypeCaller ####
#### Usage: sbatch 3_2_index_ref_genome.sh <path/to/reference_genome_directory> <ref_genome.fna/fasta>
#### Author: Wan-Ting Huang


# GATK requires a dictionary file ending in .dict and an index file ending in .fai,
# which are files that bwa index doesn't do. So to run HaplotypeCaller later, we need
# to generate .dict and .fai files.

##############################
# Arguments to pass:
workindir=reference_genome/  # Working directory. Default: 0_AclarkiiReference/
ref_genome=aclarkii_genome_chr.fna
##############################

# Check settings
echo "Working directory:	"$workindir"
Reference genome (fasta):	"$ref_genome"


"

####################
# LOAD MODULES
####################
module load openjdk/17.0.8.1_1 gatk/4.4.0.0 samtools/1.17


####################
# BEGINNNING OF SCRIPT
####################
cd $workindir

##############################
echo "Start to index "$ref_genome

# Create FASTA sequence dictionary file
gatk CreateSequenceDictionary -R $ref_genome


# Create the fasta index file
samtools faidx $ref_genome

####################
# END
####################

echo "Indexing of reference genome for GATK finished. You can now run 3_3_gatk_haplotype_call.sh"


