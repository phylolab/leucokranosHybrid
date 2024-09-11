#!/usr/bin/bash -l

#SBATCH --job-name=BWA
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --output=outscripts/BWA_%a.out
#SBATCH --error=outscripts/BWA_%a.err
#SBATCH --time=20:00:00

#### The script is used to map the reads to the reference genome ####


##############################
# Argument to pass: 
SampleList=${1-AllSamples.txt}			# Name of the file with sample names. Default: AllSamples.txt
workindir=${2-.}					# Working directory. Default: Current Directory
inputfileDir=${3-1_TrimmedReads}			# Name of the directory with the input fastq.gz files (should be in workindir). Default: 1_TrimmedReads
ReferenceGenome=${4-0_AclarkiiReference/AclarkiiGenome.Chr}	# Name (and location) of the index of the reference genome	
removeTEMP=${5-RemoveSAM}				# Do remove SAM and unsorted bam files, by default. If something else, do not remove them!

##############################


############################################################
## Relaunch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $SampleList | awk '{print $1}')


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $SampleList $workindir $inputfileDir $ReferenceGenome $removeTEMP
fi

############################################################
# Check settings
echo "Using following arguments: 

List of Samples:	"$SampleList"
Working directory:	"$workindir"
Trimmed data directory (if relative path, from "$workindir") :	"$inputfileDir"
Reference genome (index):	"$ReferenceGenome"
Remove or not 'temporary' sam and unsorted bam files:	"$removeTEMP

echo "
If settings are wrong, please look at instructions and arguments to pass

"
####################
# LOAD MODULES
####################

module load bwa/0.7.17 python/3.11.6 samtools/1.17 picard/3.0.0 bamtools/2.5.2 r/4.3.2

####################
# BEGINNNING OF SCRIPT
####################


# Get name of samples for array
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SampleList)

# Print the name of the sample
echo "Working on sample:	"$sample

# Change to working directory
cd $workindir

####################
# Run BWA
####################

mkdir -p 2B_MappingOutput
mkdir -p 2B_MappingStatistics

bwa mem -M -t 2 $ReferenceGenome \
$inputfileDir"/"$sample"_R1.Trimmed.Paired.fastq.gz" \
$inputfileDir"/"$sample"_R2.Trimmed.Paired.fastq.gz" \
> "./2B_MappingOutput/"$sample".BWA.Aclarkii.sam"


####################
# Convert to BAM
####################

samtools view -Sb -@ 2 \
"./2B_MappingOutput/"$sample".BWA.Aclarkii.sam" \
> "./2B_MappingOutput/"$sample".BWA.Aclarkii.bam"

####################
# Sort BAM (by coordinate)
####################
samtools sort -@ 2 \
"./2B_MappingOutput/"$sample".BWA.Aclarkii.bam" \
> "./2B_MappingOutput/"$sample".BWA.Aclarkii.Sorted.bam"

####################
# Index BAM
####################
samtools index "./2B_MappingOutput/"$sample".BWA.Aclarkii.Sorted.bam"



####################
# Generate additional mapping statistics
####################

bamtools stats -in "./2B_MappingOutput/"$sample".BWA.Aclarkii.Sorted.bam" \
 > "./2B_MappingStatistics/"$sample"MappingStats.txt"

picard CollectInsertSizeMetrics \
I="./2B_MappingOutput/"$sample".BWA.Aclarkii.Sorted.bam" \
O="./2B_MappingStatistics/"$sample".InsertMetrics.txt" \
H="./2B_MappingStatistics/"$sample".InsertPlot.pdf"

####################
# Remove sam and unsorted bam 
####################

if [ $removeTEMP == "RemoveSAM" ]
then
  echo "Removing sam and unsorted bam files."
  rm "./2B_MappingOutput/"$sample".BWA.Aclarkii.sam"
  rm "./2B_MappingOutput/"$sample".BWA.Aclarkii.bam"
  
fi

####################
# END
####################

echo "The script for sample "$sample" has finished"










