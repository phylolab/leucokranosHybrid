#!/bin/bash -l

#SBATCH --job-name=bwa
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --output=outscripts/bwa_%a.out
#SBATCH --error=outscripts/bwa_%a.err
#SBATCH --time=20:00:00

#### The script is used to map the reads to the reference genome ####
#### Author: Wan-Ting Huang



##############################
# Arguments to pass
sample_list=${1-all_samples.txt}			# Name of the file with sample names. Default: all_samples.txt
workindir=${2-.}					 		# Working directory. Default: current directory
input_dir=${3-1_trimmed_reads}				# Name of the directory with the input fastq.gz files (should be in workindir). Default: 1_trimmed_reads
ref_genome=${4-0_reference_genome/aclarkii_genome_chr}	# Name (and location) of the index of the reference genome	
remove_tmp=${5-remove_sam}					# Remove SAM and unsorted bam files, by default. 

##############################


############################################################
## Relaunch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $sample_list | awk '{print $1}')


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $sample_list $workindir $input_dir $ref_genome $remove_tmp
fi

############################################################
# Check settings
echo "Using following arguments: 

List of samples:	"$sample_list"
Working directory:	"$workindir"
Trimmed data directory (if relative path, from "$workindir") :	"$input_dir"
Reference genome (index):	"$ref_genome"
Remove or not 'temporary' sam and unsorted bam files:	"$remove_tmp

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
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $sample_list)

# Print the name of the sample
echo "Working on sample:	"$sample

# Change to working directory
cd $workindir

####################
# Run BWA
####################

mkdir -p 2B_MappingOutput
mkdir -p 2B_MappingStatistics

bwa mem -M -t 2 $ref_genome \
$input_dir"/"$sample"_R1.Trimmed.Paired.fastq.gz" \
$input_dir"/"$sample"_R2.Trimmed.Paired.fastq.gz" \
> "./2_bwa_mapping_out/"$sample".BWA.Aclarkii.sam"


####################
# Convert to BAM
####################

samtools view -Sb -@ 2 \
"./2_bwa_mapping_out/"$sample".BWA.Aclarkii.sam" \
> "./2_bwa_mapping_out/"$sample".BWA.Aclarkii.bam"

####################
# Sort BAM (by coordinate)
####################
samtools sort -@ 2 \
"./2_bwa_mapping_out/"$sample".BWA.Aclarkii.bam" \
> "./2_bwa_mapping_out/"$sample".BWA.Aclarkii.Sorted.bam"

####################
# Index BAM
####################
samtools index "./2_bwa_mapping_out/"$sample".BWA.Aclarkii.Sorted.bam"


####################
# Generate additional mapping statistics
####################

bamtools stats -in "./2_bwa_mapping_out/"$sample".BWA.Aclarkii.Sorted.bam" \
 > "./2_bwa_mapping_stat/"$sample"MappingStats.txt"

picard CollectInsertSizeMetrics \
I="./2_bwa_mapping_out/"$sample".BWA.Aclarkii.Sorted.bam" \
O="./2_bwa_mapping_stat/"$sample".InsertMetrics.txt" \
H="./2_bwa_mapping_stat/"$sample".InsertPlot.pdf"

####################
# Remove sam and unsorted bam 
####################

if [ $remove_tmp == "remove_sam" ]
then
  echo "Removing sam and unsorted bam files."
  rm "./2_bwa_mapping_out/"$sample".BWA.Aclarkii.sam"
  rm "./2_bwa_mapping_out/"$sample".BWA.Aclarkii.bam"
  
fi

####################
# END
####################

echo "The script for sample "$sample" has finished! You can now preprocess the BAM file before SNPs calling using 3_1_run_gatk_preprocessing.py."










