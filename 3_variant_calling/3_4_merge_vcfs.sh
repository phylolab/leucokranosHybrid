#!/bin/bash -l

#SBATCH --job-name=merge_vcfs
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --output=outscripts/merge_vcfs%a.out
#SBATCH --error=outscripts/merge_vcfs%a.err
#SBATCH --time=01:00:00

#### The script is used to merge gVCF files of all chromosomes into one gVCF file per sample ####
#### Example usage: sbatch 3_4_merge_vcfs.sh all_samples.txt . reference_genome/aclarkii_genome.fna 3_3_gatk_haplotype 3_4_calling_output .g.vcf.gz

###########################
# Arguments to pass
sample_list=${1-all_samples.txt}		# Name of the file with sample names. Default: all_samples.txt
workindir=${2-.}						# Working directory. Default: Current Directory
ref_genome=${3-reference_genome/aclarkii_genome.fna}
inputfile_dir=${4-3_3_gatk_haplotype}	# Name of the directory with the input mapping files (should be in workindir).
outputfile_dir=${5-3_4_calling_output}	# Name of the directory with the merged gVCF files (should be in workindir).
outputsuffix=${6-.g.vcf.gz}				# Suffix of output merged VCF files


############################################################

## Relaunch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $sample_list | awk '{print $1}')


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     
     exec sbatch --array=1-$nline $0 $sample_list $workindir $ref_genome $inputfile_dir $outputfile_dir $outputsuffix

fi

############################################################
# Get name of samples
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $sample_list)


############################################################


# Check settings
echo "Using following arguments: 

List of samples:	"$sample_list"

Working directory:	"$workindir"
VCF Files directory (if relative path, from "$workindir") :	"$inputfile_dir"/"$sample"

Output File: 	3_4_calling_output/"$sample".merged"$outputsuffix

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
mkdir -p $outputfile_dir


# Print the name of the sample
echo "Working on "$sample

# Change to the directory where the input files are
cd $inputfile_dir"/"$sample

# Generate a .list file listing all input files for each sample
if [ ! -s $sample".list" ]
	then
		ls | grep -v "tbi" | grep "vcf" > $sample".list"
		echo "List is created"
fi

# Combine variant calls of chromosomes
if [ ! -s "../../"$outputfile_dir"/"$sample".merged"$outputsuffix".tbi" ]
	then
		picard MergeVcfs \
		-I $sample".list" \
		-O "../../"$outputfile_dir"/"$sample".merged"$outputsuffix
		
		echo "Merging is completed!"
fi

####################
# END
####################

echo "The script for sample "$sample" has finished. Please run 3_5_genomics_db_import.sh"








