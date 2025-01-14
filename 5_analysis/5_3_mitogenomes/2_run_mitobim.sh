#!/bin/bash -l

#SBATCH -J mitobim_percula
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu 50G
#SBATCH --time 12:00:00
#SBATCH -o outscripts/%x.%a.out
#SBATCH -e outscripts/%x.%a.err
#SBATCH --export NONE

#### The script is used for baiting and iterative mapping of trimmed reads to reconstruct mitogenomes with the MITObim.pl script ####
#### Input files are the merged read files with the following format {sample_id}.extendedFrags.fastq.gz, a csv file with two columns (1-sample name, 2-associated reference mitogenome), 
#### the associated reference mitogenome in fasta format and the directory containing the ".maf" output assembly generated with the mira_assembly.sh script (should have as suffix _mira_assembly_assembly and keep the MIRA output structure)
#### Usage: sbatch --array=1-$(wc -l < all_samples.txt) run_mitobim_percula.sh <sample_list> <working/dir> <input/reads/dir> <kmer> <path/to/mira> <path/to/MITObim.pl> <path/to/output/directory> <path/to/ref/coi> <id_ref_coi>
#### Author: Sarah Schmid


##############################
# Arguments to pass
sample_list=${1-all_samples.txt}				# name of the file with sample names
workindir=${2-/scratch/sschmi13/4_mitogenomes}					 			# working directory
input_reads_dir=${3-merged_reads_flash}	# name of the directory with the input .extendedFrags.fastq.gz files (should be in workindir)
kmer=${4-31} 									#kmer value to run MITObim (default is 31 in MITObim)
mira_path=${5-./mira_4.0.2/bin/} 				# path to MIRA exe
mitobim=${6-./MITObim/MITObim.pl}				# path to MITObim.pl script
outputdir=${7-mitobim_outputs_percula}		# path to output directory
ref_coi=${8-references/percula_mtdna.fasta} #path to reference mitogenome
ref_ID=${9-percula_mtdna}

############################################################
# Check settings

echo "Using following arguments: 

Working directory: "$workindir"
File with list of all samples to process: "$sample_list"
Path to input directory with merged reads (if relative path, from "$workindir"): "$input_reads_dir"
MITObim will run with a kmer value of "$kmer"
MIRA exe (if relative path, from "$workindir"): "$mira_exe"
MITObim script (if relative path, from "$workindir"): "$mitobim"
Reference COI sequence used: "$ref_coi"
Output directory (if relative path, from "$workindir"): "$outputdir

echo "
If settings are wrong, please look at instructions and arguments to pass

"

####################
# LOAD MODULES
####################

# no module to load, but the MIRA executable PATH as well as the MITObim.pl script path should be specified in the command line

####################
# BEGINNNING OF SCRIPT
####################

# Create output directory
mkdir $workindir"/"$outputdir

# Get name of samples for array
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $workindir"/"$sample_list)

# Print the name of the sample
echo "Working on sample:	"$sample

####################
# RUNNING MITObim
####################

# Create directory for the sample
mkdir $workindir"/"$outputdir"/"$sample

# Move to the sample directory
cd $workindir"/"$outputdir"/"$sample

# Run MITObim command
$workindir"/"$mitobim -start 1 -end 60 \
	-sample $sample \
	-ref $ref_ID \
	-readpool $workindir"/"$input_reads_dir"/"$sample.extendedFrags.fastq.gz \
	--quick $workindir"/"$ref_coi \
	--clean \
	--kbait $kmer \
	--mirapath $workindir"/"$mira_path

####################
# END
####################

echo "The baiting and iterative mapping for sample "$sample" is done!"

