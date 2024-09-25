#!/usr/bin/bash -l

#SBATCH --job-name=ELAI_vcf2bimbam
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=outscripts/elai_vcf2bimbam.out
#SBATCH --error=outscripts/elai_vcf2bimbam.err
#SBATCH --time=1:00:00

#### The script is used to convert VCF file to bimbam file which is the required file format for ELAI analysis ####

###########################
# Arguments to pass:
workindir=5_ELAI			# Working directory. Default: 5_ELAI
inputfileDir=0_annotatedvcf		# Name of the directory with annotated files. Default: 0_annotatedvcf
outputfileDir=1_bimbam			# Name of the directory with output files. Default: 1_bimbam
poplistDir=pop_files			# Name of the directory with files listing out the samples of each population. Default: pop_files

############################################################


####################
# LOAD MODULES
####################

module load plink/1.9-beta6.27


####################
# BEGINNING OF SCRIPT
####################
cd $workindir

mkdir -p 1_bimbam

for file in $inputfileDir/*.annotated.vcf
do
	echo "Working on "$file
	
	string=$(cut -d/ -f2 <<< $file)
	chr=$(cut -d_ -f1 <<< $string)
	
	for pop in `cat $poplistDir/pop.list`
	do
		plink --vcf $file \
			--recode bimbam --keep $poplistDir/$pop \
			--out $outputfileDir"/"$chr"_SNP.SG.Filtered.nospandel.${pop%%.pop}"
	done
	
	echo "Done for "$file"
	-"
	
done
		
echo "
--------------------

The script has finished."	
		
		
		
