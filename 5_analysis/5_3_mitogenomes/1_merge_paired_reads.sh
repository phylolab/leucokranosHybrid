#!/bin/sh

#SBATCH -J merge_reads
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu 10G
#SBATCH --time 20:00:00
#SBATCH -o outscripts/%x.%a.%A.out
#SBATCH -e outscripts/%x.%a.%A.err
#SBATCH --export NONE


#### The script is used to merge overlapping paired end reads with FLASh ####
#### Inputs are trimmed paired reads with the following format: [sample_id]_R[1-2].Trimmed.Paired.fastq.gz
#### Usage with slurm: sbatch --array=1-$(wc -l < <sample_id_list>.txt) merge_paired_reads.sh <sample_id_list.txt> <working/directory> <input/directory/with/fastq> <path/to/output/directory> <path/to/flash/exe>
#### Author: Sarah Schmid

##############################
# Arguments to pass

sample_list=${1-all_samples.txt} # name of the file with sample names (e.g. 388_WN). 
workindir=${2-.} # working directory. 
input_dir=${3-fastq_files/} # name of the directory with the input fastq files (should be in workindir).
output_dir=${4-merged_reads_flash/} # name of the output directory
flash=${5-./flash_1.2.11/flash} # FLASh executable

############################################################
# Check settings
echo "Using following arguments: 

List of samples:	"$sample_list"
Working directory:	"$workindir"
Trimmed data directory (if relative path, from "$workindir") :	"$input_dir"
FLASh executable path: "$flash

echo "
If settings are wrong, please look at instructions and arguments to pass

"

####################
# LOAD MODULES
####################

# no module to load, but the FLASh executable PATH should be specified in the command line

####################
# BEGINNNING OF SCRIPT
####################

# Change to working directory and create output directory
cd $workindir
mkdir -p $output_dir


# Get name of samples for array
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $workindir"/"$sample_list)

# Print the name of the sample
echo "Working on sample:	"$sample


####################
# RUNNING FLASH
####################

cd $workindir

$flash $input_dir"/"$sample"_R1.Trimmed.Paired.fastq"  $input_dir"/"$sample"_R2.Trimmed.Paired.fastq" --compress --max-overlap=120 --allow-outies --interleaved-output --output-prefix=$sample --output-directory=$output_dir &> "flash_log/"$sample"_flash.log"

####################
# END
####################

echo "The script for sample "$sample" has finished! You can now run the script mitobim_full_mito.sh to reconstruct the mitogenome. "



