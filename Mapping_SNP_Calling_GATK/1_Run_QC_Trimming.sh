#!/usr/bin/bash -l

#SBATCH --job-name=QCTrim
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --output=outscripts/QCTrim_%a.out
#SBATCH --error=outscripts/QCTrim_%a.err
#SBATCH --time=16:00:00

#### The script is used to perform quality tests and trimming on sequencing data ####

##############################
# Argument to pass: 
SampleList=${1-AllSamples.txt}		# Name of the file with sample names. (should be in workindir) Default: AllSamples.txt
workindir=${2-.}				# Working directory. Default: Current Directory
inputfileDir=${3-/Path/To/RawData}			# Name of the directory with the input fastq.gz files (should be in workindir). Default: /Path/To/RawData
						# The fastq files should be in the format $sample_R1.fastq.gz and $sample_R2.fastq.gz
adapterfiles=${4-Mapping_SNP_Calling_GATK/TruSeq3-PE-2.fa}	# Path and name of the files where adapter sequences are. Default: Mapping_SNP_Calling_GATK/TruSeq3-PE-2.fa
minQuality=${5-3}				# Minimum quality score for trimming. Default: 20
minLength=${6-36}				# Minimum length of reads: Default: 50
headcrop=${7-0}
##############################

##############################
## Relaunch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $SampleList | awk '{print $1}')


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $SampleList $workindir $inputfileDir $adapterfiles $minQuality $minLength $headcrop
fi

##############################
# Check settings
echo "Using following arguments: 

List of Samples:	"$SampleList"
Working directory:	"$workindir"
Raw data directory (if relative path, from "$workindir" :	"$inputfileDir"
Adapter files: 	"$adapterfiles"
Minimum quality of bases for trimming:	"$minQuality"
Minimum length of reads:	"$minLength"
Number of basis cut beginning of reads:	"$headcrop

echo "
If settings are wrong, please look at instructions and arguments to pass
"

module load fastqc/0.12.1 trimmomatic/0.39

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
# Run quality report with fastqc for each read file 
####################

# make directory fore fastqc raw data if does not exist
mkdir -p 1_fastQC_RawData

# generate fastqc report
fastqc -o 1_fastQC_RawData $inputfileDir"/"$sample"_R1.fastq.gz"
fastqc -o 1_fastQC_RawData $inputfileDir"/"$sample"_R2.fastq.gz"

####################
# Run the trimming and adapter check with trimmomatic
####################

# make directory fore Trimmed reads and Summary of trimming if does not exist
mkdir -p 1_TrimmedReads
mkdir -p 1_TrimmingStatistics

trimmomatic PE -phred33 -summary "1_TrimmingStatistics/"$sample".TrimmingStats.txt" \
$inputfileDir"/"$sample"_R1.fastq.gz" $inputfileDir"/"$sample"_R2.fastq.gz" \
"1_TrimmedReads/"$sample"_R1.Trimmed.Paired.fastq.gz" "1_TrimmedReads/"$sample"_R1.Trimmed.Unpaired.fastq.gz" \
"1_TrimmedReads/"$sample"_R2.Trimmed.Paired.fastq.gz" "1_TrimmedReads/"$sample"_R2.Trimmed.Unpaired.fastq.gz" \
ILLUMINACLIP:$adapterfiles:2:30:10:2:keepBothReads SLIDINGWINDOW:4:15 \
LEADING:$minQuality TRAILING:$minQuality MINLEN:$minLength HEADCROP:$headcrop


####################
# Run quality report with fastqc for each Trimmed read file 
####################
mkdir -p fastQC_TrimmedData
fastqc -o fastQC_TrimmedData "1_TrimmedReads/"$sample"_R1.Trimmed.Paired.fastq.gz"
fastqc -o fastQC_TrimmedData "1_TrimmedReads/"$sample"_R2.Trimmed.Paired.fastq.gz"
fastqc -o fastQC_TrimmedData "1_TrimmedReads/"$sample"_R1.Trimmed.Unpaired.fastq.gz"
fastqc -o fastQC_TrimmedData "1_TrimmedReads/"$sample"_R2.Trimmed.Unpaired.fastq.gz"


####################
# END
####################

echo "The script for sample "$sample" has finished"




