#!/bin/bash -l

#SBATCH --job-name=gatk_haplo_call
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --output=outscripts/gatk_haplo_call_1_%a.out
#SBATCH --error=outscripts/gatk_haplo_call_1_%a.err
#SBATCH --time=72:00:00

#### The script is used to call variants for each sample by chromosome using GATK HaplotypeCaller ####
#### Example usage: sbatch 3_3_gatk_haplotype_call.sh all_samples.txt . reference_genome/aclarkii_genome_chr.fna 3_1_gatk_preprocessing .dedupped.reorder 3_3_gatk_haplotype .g.vcf.gz reference_genome/chromosomes.txt
#### Author: Wan-Ting Huang

#!! Change .out and .err name according to the number of round you are doing 

##############################
# Argument to pass
sample_list=${1-all_samples.txt}		# Name of the file with sample names. Default: all_samples.txt
workindir=${2-.}						# Working directory. Default: Current Directory
ref_genome=${3-reference_genome/aclarkii_genome_chr.fna}
inputfile_dir=${4-3_1_gatk_preprocessing}	# Name of the directory with the pre-processed mapping files (should be in workindir).
inputsuffix=${5-.dedupped.reorder}			# Suffix of input bam files
outputfile_dir=${6-3_3_gatk_haplotype}		# Name of the directory with the output gVCF files (should be in workindir).
outputsuffix=${7-.g.vcf.gz}					# Suffix of output intermediate gVCF files
file_chromosomes=${8-reference_genome/chromosomes.txt}   # File with chromosomes names
##############################


############################################################
## Relaunch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $sample_list | awk '{print $1}')


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $sample_list $workindir $ref_genome $inputfile_dir $inputsuffix $outputfile_dir $outputsuffix $file_chromosomes
fi

############################################################
# Check settings
# Get name of samples for array
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $sample_list)

echo "Using following arguments: 

List of Samples:	"$sample_list"
Working directory:	"$workindir"

Filtered Mapping Files directory (if relative path, from "$workindir") :	"$inputfile_dir"
Input File:	"$sample$inputsuffix".bam

gVCF Files directory (if relative path, from "$workindir") :	"$outputfile_dir"
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

mkdir -p $inputfile_dir"/"$sample


if [ -z "$(ls -A $inputfile_dir"/"$sample)" ] # If the folder $inputfile_dir/$sample is empty 
											#(the wrapped double quotation marks here are essential)
then										# then split the BAM file
	echo "** Split merged reads BAM file **"
	cp $inputfile_dir"/"$sample$inputsuffix".bam" $inputfile_dir"/"$sample"/"$sample$inputsuffix".bam"
	bamtools split -in $inputfile_dir"/"$sample"/"$sample$inputsuffix".bam" -reference
fi


####################
# Remove gVCF file with unfinished calling
####################
#!! Change error file name according to the number of round you are doing 
if [[ -n $(grep "TIME LIMIT" outscripts/gatk_haplo_call_${SLURM_ARRAY_TASK_ID}.err) ]]
			# If the output of grep command has a non-zero length
	then 	# then grep the name of chromosome that is unfinished
		chr=$(tail -n3 outscripts/gatk_haplo_call_${SLURM_ARRAY_TASK_ID}.err | grep "Chr.*:" -o | cut -d: -f1 | uniq)
		if [[ -n ${chr} ]] # Test whether we grep the chromosome name
			then
				echo "!!"$sample"."$chr$outputsuffix" did not finish properly due to time limit."
				rm $outputfile_dir"/"$sample"/"$sample"."$chr$outputsuffix
				
				echo "Remove the generated file and re-call on this chromosome.
				-----"
		fi
fi

####################
# Create index files and generate gVCF files
####################
mkdir -p $outputfile_dir"/"$sample

while read -r chrid; do

if [ ! -s $inputfile_dir"/"$sample"/"$sample$inputsuffix".REF_"$chrid".bam.bai" ] 
	# If there is no this file or the file is empty (-s test for an existence of non-empty files)
	then 		# then index the BAM file
		echo "** Index "$chrid" BAM file **"
		samtools index $inputfile_dir"/"$sample"/"$sample$inputsuffix".REF_"$chrid".bam"
fi

if [ ! -s $outputfile_dir"/"$sample"/"$sample"."$chrid$outputsuffix".tbi" ]
	# If the gVCF file for this chromosome does not exist
	then	# then call on this chromosome
		
		echo "...Call on "$chrid
		
		gatk --java-options "-Xmx4g" HaplotypeCaller -R $ref_genome \
		-I $inputfile_dir"/"$sample"/"$sample$inputsuffix".REF_"$chrid".bam" \
		-O $outputfile_dir"/"$sample"/"$sample"."$chrid$outputsuffix \
		-L $chrid \
		-ERC GVCF
fi

done < $file_chromosomes

####################
# END
####################

echo "The script for sample "$sample" has finished, you can now run 3_4_merge_vcfs.sh"
