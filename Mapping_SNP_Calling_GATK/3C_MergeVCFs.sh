#!/usr/bin/bash -l

#SBATCH --job-name=MergeVcfs
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --output=outscripts/MergeVcfs_%a.out
#SBATCH --error=outscripts/MergeVcfs_%a.err
#SBATCH --time=01:00:00

#### The script is used to merge gVCF files of all chromosomes into one gVCF file per sample ####


###########################
# Arguments to pass:
SampleList=${1-AllSamples.txt}	# Name of the file with sample names. Default: AllSamples.txt
workindir=${2-.}					# Working directory. Default: Current Directory
ReferenceGenome=${3-0_AclarkiiReference/AclarkiiGenome.Chr.fna}
inputfileDir=${4-3B_CallingHaplotype}			# Name of the directory with the input mapping files (should be in workindir). Default: 3B_CallingHaplotype
outputfileDir=${5-3C_CallingOutput}			# Name of the directory with the merged gVCF files (should be in workindir). Default: 3C_CallingOutput
outputsuffix=${6-.g.vcf.gz}	# Suffix of output merged VCF files



############################################################

## Relaunch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $SampleList | awk '{print $1}')


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     
     exec sbatch --array=1-$nline $0 $SampleList $workindir $ReferenceGenome $inputfileDir $outputfileDir $outputsuffix

fi

############################################################
# Get name of samples
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SampleList)


############################################################


# Check settings
echo "Using following arguments: 

List of Samples:	"$SampleList"

Working directory:	"$workindir"
VCF Files directory (if relative path, from "$workindir") :	"$inputfileDir"/"$sample"

Output File: 	"$outputfileDir"/"$sample".Merged"$outputsuffix

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
mkdir -p 3C_CallingOutput


# Print the name of the sample
echo "Working on "$sample

# Change to the directory where the input files are
cd $inputfileDir"/"$sample

# Generate a .list file listing all input files for each sample
if [ ! -s $sample".list" ]
	then
		ls | grep -v "tbi" | grep "vcf" > $sample".list"
		echo "List is created"
fi

# Combine variant calls of chromosomes
if [ ! -s "../../"$outputfileDir"/"$sample".Merged"$outputsuffix".tbi" ]
	then
		picard MergeVcfs \
		-I $sample".list" \
		-O "../../"$outputfileDir"/"$sample".Merged"$outputsuffix
		
		echo "Merging is completed!"
fi

####################
# END
####################

echo "The script for sample "$sample" has finished"








