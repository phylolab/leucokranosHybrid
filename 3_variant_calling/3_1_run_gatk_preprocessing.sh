#!/bin/bash -l

#SBATCH --job-name=gatk_preprocess
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=10G
#SBATCH --output=outscripts/gatk_preproc_%a.out
#SBATCH --error=outscripts/gatk_preproc_%a.err
#SBATCH --time=20:00:00

#### The script is used to process mapping files before variant calling with GATK ####
#### Companion script of 3_1_gatk_preprocessing.sh
#### Input bam file names should have the following format SAMPLEID.BWA.Aclarkii.Sorted.bam 
#### Command line example of how to run this script in slurm: sbatch 3_1_run_gatk_preprocessing.sh all_samples.txt /scratch/sschmi13/leucokranos/2_bwa_mapping/ 2_bwa_mapping_out/ reference_genome/aclarkii_genome_chr.fna 3 102400000
#### Author: Sarah Schmid


##############################
# Arguments to pass
sample_list=${1-all_samples.txt}			# Name of the file with sample names. Default: all_samples.txt
workindir=${2-.}					 		# Working directory. Default: current directory (should end with a "/")
input_dir=${3-2_bwa_mapping_out}			# Name of the directory with output bam files (should be in workindir). Default: 2_bwa_mapping_out
ref_genome=${4-reference_genome/aclarkii_genome_chr}	# Name (and location) of the index of the reference genome	
max_cpu=${5-3}								# Maximum number of CPU to use. Default = 3
max_mem=${6-102400000}						# Maximum memory to use. Default = 102400000

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
Max computer:	"$max_cpu"
Max memory: "$max_mem

echo "
If settings are wrong, please look at instructions and arguments to pass

"
####################
# LOAD MODULES
####################

module load python/3.11.6 

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
# Running GATK preprocessing python script
####################

python 3_1_gatk_preprocessing.py $workindir $input_dir"/"$sample".BWA.Aclarkii.Sorted.bam" $sample $max_cpu $max_mem $ref_genome

####################
# END
####################

echo "The script for sample "$sample" has finished! You can now either do base recalibration with 3_2_gatk_base_recalibration.py or directly proceed to snp calling with 3_3_gatk_snp_calling.py "



