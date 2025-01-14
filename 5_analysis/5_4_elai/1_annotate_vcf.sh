#!/usr/bin/bash -l

#SBATCH --job-name=elai_annotate
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=outscripts/%x.%a.out
#SBATCH --error=outscripts/%x.%a.err
#SBATCH --time=1:00:00

### The script is used to prepare annotated VCF file for subsequent conversion into bimbam file needed in ELAI ###
### Author: Wan-Ting Huang


###########################
# Arguments to pass:
workindir=${1-/scratch/sschmi13/5_elai}				
inputdir=${2-/scratch/sschmi13/final_vcf/chryso_san_leuco_vcfs/} # directory with VCF files 
outputdir=${3-annotated_vcf}			
chrlist=${4-chr.list} # chromosome list matching format in VCF file

############################################################


####################
# LOAD MODULES
####################

module load plink/1.9-beta6.27 bcftools/1.16


####################
# BEGINNING OF SCRIPT
####################
cd $workindir

mkdir -p $outputdir

while read -r chr
do
	bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
			$inputdir"/"$chr".no_lowdp_indv_chryso_sand_leuco.gatk_filtered.no_mono.biall.mac.maxmiss.mindp.maxdp.qual.vcf.gz" \
			> $outputdir"/"$chr".no_lowdp_indv_chryso_sand_leuco.gatk_filtered.no_mono.biall.mac.maxmiss.mindp.maxdp.qual.annotated.vcf"
			
	echo "Done for "$chr
	
done < $chrlist