#!/bin/bash -l

#SBATCH --job-name=qc_trim
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --output=outscripts/qc_trim_%a.out
#SBATCH --error=outscripts/qc_trim_%a.err
#SBATCH --time=16:00:00

#### The script is used to perform quality tests and trimming of raw sequencing reads ####

##############################
## Arguments to pass 
samples_list=${1-all_samples.txt}		# Name of the file with sample names. (should be in workindir) Default: all_samples.txt
workindir=${2-.}						# Working directory. Default: Current Directory
inputfiledir=${3-/path/to/rawdata}		# Name of the directory with the input $sample_R1/2_fastq.gz files (should be in workindir). Default: /path/to/rawdata								
adapterfiles=${4-TruSeq3-PE-2.fa}	    # Name of the files where adapter sequences are. Default: TruSeq3-PE-2.fa
min_quality=${5-3}						# Minimum quality score for trimming. Default: 20
min_length=${6-36}						# Minimum length of reads: Default: 50
headcrop=${7-0}
##############################

##############################
## Relaunch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $samples_list | awk '{print $1}')


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $samples_list $workindir $inputfiledir $adapterfiles $min_quality $min_length $headcrop
fi

##############################
# Check settings
echo "Using following arguments: 

List of Samples:	"$samples_list"
Working directory:	"$workindir"
Raw data directory (if relative path, from "$workindir" :	"$inputfiledir"
Adapter files: 	"$adapterfiles"
Minimum quality of bases for trimming:	"$min_quality"
Minimum length of reads:	"$min_length"
Number of basis cut beginning of reads:	"$headcrop

echo "
If settings are wrong, please look at instructions and arguments to pass
"

module load fastqc/0.12.1 trimmomatic/0.39

####################
# BEGINNNING OF SCRIPT
####################

# Get name of samples for array
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $samples_list)

# Print the name of the sample
echo "Working on sample:	"$sample

# Change to working directory
cd $workindir

####################
# Run quality report with fastqc for each read file 
####################

# make directory fore fastqc raw data if does not exist
mkdir -p 1_fastqc_rawdata

# generate fastqc report
fastqc -o 1_fastqc_rawdata $inputfiledir"/"$sample"_R1.fastq.gz"
fastqc -o 1_fastqc_rawdata $inputfiledir"/"$sample"_R2.fastq.gz"

####################
# Run the trimming and adapter check with trimmomatic
####################

# make directory fore Trimmed reads and Summary of trimming if does not exist
mkdir -p 1_trimmed_reads
mkdir -p 1_trimming_statistics

trimmomatic PE -phred33 -summary "1_trimming_statistics/"$sample".trimming_stats.txt" \
$inputfiledir"/"$sample"_R1.fastq.gz" $inputfiledir"/"$sample"_R2.fastq.gz" \
"1_trimmed_reads/"$sample"_R1.trimmed.paired.fastq.gz" "1_trimmed_reads/"$sample"_R1.trimmed.unpaired.fastq.gz" \
"1_trimmed_reads/"$sample"_R2.trimmed.paired.fastq.gz" "1_trimmed_reads/"$sample"_R2.trimmed.unpaired.fastq.gz" \
ILLUMINACLIP:$adapterfiles:2:30:10:2:keepBothReads SLIDINGWINDOW:4:15 \
LEADING:$min_quality TRAILING:$min_quality MINLEN:$min_length HEADCROP:$headcrop


####################
# Run quality report with fastqc for each trimmed read file 
####################
mkdir -p 1_fastqc_trimmed_data
fastqc -o 1_fastqc_trimmed_data "1_trimmed_reads/"$sample"_R1.trimmed.paired.fastq.gz"
fastqc -o 1_fastqc_trimmed_data "1_trimmed_reads/"$sample"_R2.trimmed.paired.fastq.gz"
fastqc -o 1_fastqc_trimmed_data "1_trimmed_reads/"$sample"_R1.trimmed.unpaired.fastq.gz"
fastqc -o 1_fastqc_trimmed_data "1_trimmed_reads/"$sample"_R2.trimmed.unpaired.fastq.gz"


####################
# END
####################

echo "The script for sample "$sample" has finished"



