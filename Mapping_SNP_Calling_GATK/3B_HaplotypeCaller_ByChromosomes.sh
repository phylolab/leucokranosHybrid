#!/usr/bin/bash -l

#SBATCH --job-name=HaplotypeCaller
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --output=outscripts/Haplotype1_%a.out
#SBATCH --error=outscripts/Haplotype1_%a.err
#SBATCH --time=72:00:00

#### The script is used to call variants for each sample by chromosome using GATK HaplotypeCaller ####

#!! Change .out and .err name according to the number of round you are doing 

##############################
# Argument to pass: 
SampleList=${1-AllSamples.txt}		# Name of the file with sample names. Default: AllSamples.txt
workindir=${2-.}				# Working directory. Default: Current Directory
ReferenceGenome=${3-0_AclarkiiReference/AclarkiiGenome.Chr.fna}
inputfileDir=${4-2C_FilteredMapping}			# Name of the directory with the input mapping files (should be in workindir). Default: 2C_FilteredMapping
inputsuffix=${5-.BWA.Aclarkii.Sort.Filt_mergedReads}		# Suffix of input bam files
outputfileDir=${6-3B_CallingHaplotype}			# Name of the directory with the output gVCF files (should be in workindir). Default: 3B_CallingHaplotype
outputsuffix=${7-.g.vcf.gz}	# Suffix of output intermediate gVCF files
FileChromosomes=${8-0_AclarkiiReference/Chromosomes.txt}   # File with chromosomes names
##############################


############################################################
## Relounch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $SampleList | awk '{print $1}')


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $SampleList $workindir $ReferenceGenome $inputfileDir $inputsuffix $outputfileDir $outputsuffix $FileChromosomes
fi

############################################################
# Check settings
# Get name of samples for array
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SampleList)

echo "Using following arguments: 

List of Samples:	"$SampleList"
Working directory:	"$workindir"

Filtered Mapping Files directory (if relative path, from "$workindir") :	"$inputfileDir"
Input File:	"$sample$inputsuffix".bam

gVCF Files directory (if relative path, from "$workindir") :	"$outputfileDir"
Output File: 	"$sample$outputsuffix

echo "
If settings are wrong, please look at instructions and arguments to pass

"


####################
# LOAD MODULES
####################

module load openjdk/17.0.8.1_1 gatk/4.4.0.0 samtools/1.17 bamtools/2.5.2

####################
# BEGINNNING OF SCRIPT
####################
# Print the name of the sample
echo "Working on sample:	"$sample


# Change to working directory
cd $workindir

####################
# SPLIT BAM by Chromosome
####################

mkdir -p $inputfileDir"/"$sample


if [ -z "$(ls -A $inputfileDir"/"$sample)" ] # If the folder $inputfileDir/$sample is empty 
											#(the wrapped double quotation marks here are essential)
then										# then split the BAM file
	echo "** Split merged reads BAM file **"
	cp $inputfileDir"/"$sample$inputsuffix".bam" $inputfileDir"/"$sample"/"$sample$inputsuffix".bam"
	bamtools split -in $inputfileDir"/"$sample"/"$sample$inputsuffix".bam" -reference
fi


####################
# Remove gVCF file with unfinished calling
####################
#!! Change error file name according to the number of round you are doing 
if [[ -n $(grep "TIME LIMIT" outscripts/Haplotype_${SLURM_ARRAY_TASK_ID}.err) ]]
			# If the output of grep command has a non-zero length
	then 	# then grep the name of chromosome that is unfinished
		chr=$(tail -n3 outscripts/Haplotype_${SLURM_ARRAY_TASK_ID}.err | grep "Chr.*:" -o | cut -d: -f1 | uniq)
		if [[ -n ${chr} ]] # Test whether we grep the chromosome name
			then
				echo "!!"$sample"."$chr$outputsuffix" did not finish properly due to time limit."
				rm $outputfileDir"/"$sample"/"$sample"."$chr$outputsuffix
				
				echo "Remove the generated file and re-call on this chromosome.
				-----"
		fi
fi

####################
# Create index files and Generate GVCF files
####################
mkdir -p $outputfileDir"/"$sample

while read -r chrid; do

if [ ! -s $inputfileDir"/"$sample"/"$sample$inputsuffix".REF_"$chrid".bam.bai" ] 
	# If there is no this file or the file is empty (-s test for an existence of non-empty files)
	then 		# then index the BAM file
		echo "** Index "$chrid" BAM file **"
		samtools index $inputfileDir"/"$sample"/"$sample$inputsuffix".REF_"$chrid".bam"
fi

if [ ! -s $outputfileDir"/"$sample"/"$sample"."$chrid$outputsuffix ]
	# If the gVCF file for this chromosome does not exist
	then	# then call on this chromosome
		
		echo "...Call on "$chrid
		
		gatk --java-options "-Xmx4g" HaplotypeCaller -R $ReferenceGenome \
		-I $inputfileDir"/"$sample"/"$sample$inputsuffix".REF_"$chrid".bam" \
		-O $outputfileDir"/"$sample"/"$sample"."$chrid$outputsuffix \
		-L $chrid \
		-ERC GVCF
fi

done < $FileChromosomes

####################
# END
####################

echo "The script for sample "$sample" has finished"
