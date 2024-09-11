#!/usr/bin/bash -l

#SBATCH --job-name=GenomicsDBImport
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --output=outscripts/GenomicsDBImport.out
#SBATCH --error=outscripts/GenomicsDBImport.err
#SBATCH --time=40:00:00

#### The script is used to import all the gVCF files into a database workspace which would facilitate the process of joint genotyping ####

###########################
# Arguments to pass:
workindir=.					# Working directory. Default: Current directory
ReferenceGenome=0_AclarkiiReference/AclarkiiGenome.Chr.fna
inputfileDir=3C_CallingOutput			# Name of the directory with the merged gVCF files each sample (should be in workindir). Default: 3C_CallingOutput
outputfileDir=/scratch/whuang4/GenomicsDB			# Name of the directory with a database workspace. Default: /scratch/whuang4/GenomicsDB
FileChromosomes=0_AclarkiiReference/Chromosomes.list

############################################################


# Check settings
echo "Using following arguments: 

Working directory:	"$workindir"
Merged gVCF Files directory (if relative path, from "$workindir") :	"$inputfileDir"
Output database: "$outputfileDir


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

cd $workindir

# Generate a .sample_map file listing all input files
for i in ${inputfileDir}/*.vcf.gz
do
		sample_name=$(basename ${i} .Merged.g.vcf.gz)
		file_name=$(echo "$i")
		echo -e "$sample_name\t$file_name"
		
done > $inputfileDir"/AllSamples.sample_map"


echo "List is created"

# Import variant calls of multiple samples by chromosome to a database
echo "------------------
Start to import..."

gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
	--genomicsdb-workspace-path "${outputfileDir}" \	# The wrapped double quotes are essentail here
	--sample-name-map $inputfileDir"/AllSamples.sample_map" \
	-L $FileChromosomes
	
echo "Importation finished. Copy the database to work directory."

# Copy the database folder to /work
cp -r $outputfileDir $workindir


####################
# END
####################

				
echo "The script has finished!"				
	



