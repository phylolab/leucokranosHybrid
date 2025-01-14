#!/bin/bash

#SBATCH --job-name=GenomicsDBImport
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --output=outscripts/GenomicsDBImport.out
#SBATCH --error=outscripts/GenomicsDBImport.err
#SBATCH --time=40:00:00
#SBATCH --mail-user sarah.schmid@unil.ch
#SBATCH --mail-type ALL


#### The script is used to import all the gVCF files into a database workspace which would facilitate the process of joint genotyping ####
#### Author: Wan-Ting Huang


###########################
# Arguments to pass
workindir=${1-/scratch/sschmi13/leucokranos/2_bwa_mapping}					
ref_genome=${2-reference_genome/aclarkii_genome_chr.fna}
inputfile_dir=${3-3_4_calling_output}			                  	# Name of the directory with the merged gVCF files each sample (should be in workindir).
outputfile_dir=${4-/scratch/sschmi13/leucokranos/GenomicsDB}	 	# Name of the directory with a database workspace.
file_chromosomes=${5-reference_genome/chromosomes.list}   			# File with chromosomes names

############################################################


# Check settings
echo "Using following arguments: 

Working directory:	"$workindir"
Merged gVCF Files directory (if relative path, from "$workindir") :	"$inputfile_dir"
Output database: "$outputfile_dir


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

# Generate a .sample_map file listing all input files

# Remove file which is already present
rm $workindir"/"$inputfile_dir"/all_samples.sample_map"

for i in ${workindir}/${inputfile_dir}/*.vcf.gz
do
	sample_name=$(basename ${i} .merged.g.vcf.gz)
	file_name=$(echo "$i")
	echo -e "$sample_name\t$file_name" >> $workindir"/"$inputfile_dir"/all_samples.sample_map"
		
done 


echo "List is created"

# Import variant calls of multiple samples by chromosome to a database
echo "------------------
Start to import..."

gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
	--genomicsdb-workspace-path "${outputfile_dir}" \
	--sample-name-map $workindir"/"$inputfile_dir"/all_samples.sample_map" \
	-L $workindir"/"$file_chromosomes
	
echo "Importation finished. Copy the database to work directory."

# Copy the database folder to /work
cp -r $outputfile_dir $workindir


####################
# END
####################

				
echo "The script has finished! Please run 3_6_joint_genotyping.sh"				
	



