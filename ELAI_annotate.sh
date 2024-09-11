#!/usr/bin/bash -l

#SBATCH --job-name=ELAI_annotate
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=outscripts/elai_annotate.out
#SBATCH --error=outscripts/elai_annotate.err
#SBATCH --time=1:00:00

### The script is used to prepare annotated VCF file for subsequent conversion into bimbam file needed in ELAI ###

###########################
# Arguments to pass:
workindir=.					# Working directory. Default: Current directory
inputfileDir=4B_FilteredCalling			# Name of the directory with the filtered callset (should be in workindir). Default: 4B_FilteredCalling
outputfileDir=5_ELAI/0_annotatedvcf			# Default: 5_ELAI/0_annotatedvcf
chrlist=0_AclarkiiReference/Chromosomes.list

############################################################


####################
# LOAD MODULES
####################

module load plink/1.9-beta6.27 bcftools/1.16



####################
# BEGINNING OF SCRIPT
####################
cd $workindir

mkdir -p 5_ELAI/0_annotatedvcf

while read -r chr
do
	bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
			$inputfileDir"/SpeciationGenomicsTutorial/"$chr"_SNP.SG.Filtered.nospandel.vcf.gz" \
			> $outputfileDir"/"$chr"_SNP.SG.Filtered.nospandel.annotated.vcf"
			
	echo "Done for "$chr
	
done < $chrlist