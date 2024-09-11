#!/usr/bin/bash -l

#SBATCH --job-name=IndexGATK
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --output=outscripts/IndexGATK.out
#SBATCH --error=outscripts/IndexGATK.err
#SBATCH --time=01:00:00

#### The script is used to generate index files for performing GATK HaplotypeCaller ####

# GATK requires a dictionary file ending in .dict and an index file ending in .fai,
# which are files that bwa index doesn't do. So to run HaplotypeCaller later, we need
# to generate .dict and .fai files.

##############################
# Argument to pass: 
workindir=0_AclarkiiReference  # Working directory. Default: 0_AclarkiiReference/
ReferenceGenome=AclarkiiGenome.Chr.fna
##############################

# Check settings
echo "Working directory:	"$workindir"
Reference genome (fasta):	"$ReferenceGenome"


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
echo "Start to index "$ReferenceGenome

# Create FASTA sequence dictionary file
gatk CreateSequenceDictionary -R $ReferenceGenome


# Create the fasta index file
samtools faidx $ReferenceGenome

####################
# END
####################

echo "Indexing finished"


