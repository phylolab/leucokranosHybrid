#!/usr/bin/bash -l

#SBATCH --job-name=JointGenotyping
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --output=outscripts/Genotyping.out
#SBATCH --error=outscripts/Genotyping.err
#SBATCH --time=50:00:00

#### The script is used to perform GATK GenotypeGVCFs ####

###########################
# Arguments to pass:
workindir=.					# Working directory. Default: Current directory
ReferenceGenome=0_AclarkiiReference/AclarkiiGenome.Chr.fna
db=GenomicsDB			# Name of database (should be in workindir). Need to copy from /scratch to working directory since it wouldn't work with path other than gendb:// Default: GenomicsDB
outputfileDir=3C_CallingOutput	# Name of the directory with a vcf file containing genotypes of all samples. Default: 3C_CallingOutput

############################################################


# Check settings
echo "Using following arguments: 

Working directory:	"$workindir"
Name of database :	"$db"
Output directory: "$outputfileDir


echo "
If settings are wrong, please look at instructions and arguments to pass

"

####################
# LOAD MODULES
####################

module load openjdk/17.0.8.1_1 gatk/4.4.0.0 picard/3.0.0 samtools/1.17 bamtools/2.5.2 

####################
# BEGINNING OF SCRIPT
####################

cd $workindir


# Perform joint genotyping on GenomicsDB workspace 
echo "Start to joint genotyping..."

gatk --java-options "-Xmx4g" GenotypeGVCFs \
	-R $ReferenceGenome \
	-V gendb://$db \
	-O $outputfileDir"/AllSamples.vcf.gz"
	
echo "Joint genotyping is done and output a file "$outputfileDir"/AllSamples.vcf.gz."

####################
# END
####################

echo "The script has finished!"


