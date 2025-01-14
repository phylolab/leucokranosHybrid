#!/bin/bash

#SBATCH --job-name=gathervcfs
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --output=outscripts/gathervcfs.out
#SBATCH --error=outscripts/gathervcfs.err
#SBATCH --time=02:00:00

#### The script is used to perform PICARD GatherVcfs to create a single VCF file with all chromosomes and samples included. ####
#### List of chromosomes to merge hardcoded, please change directly in the script according to your project. 


############################################################
# Arguments to pass
inputdir=${1-/scratch/sschmi13/leucokranos/2_bwa_mapping/3_6_joint_geno_out} # Directory with all the individual chromosomes VCF files			
output=${2-all_chr_all_samples.vcf.gz} # Output file name


############################################################
# Check settings
echo "Using following arguments: 

Working directory:	"$inputdir"
Output file: "$output

echo "
Important: list of chromosomes to merge is hardcoded, please change directly in the script (lines 44 to 49) according to your project. 

"

####################
# LOAD MODULES
####################

module load openjdk/17.0.8.1_1 picard/3.0.0

####################
# BEGINNING OF SCRIPT
####################

picard GatherVcfs I=$inputdir"/Chr1_all_samples.vcf.gz" I=$inputdir"/Chr2_all_samples.vcf.gz" I=$inputdir"/Chr3_all_samples.vcf.gz" I=$inputdir"/Chr4_all_samples.vcf.gz" \
I=$inputdir"/Chr5_all_samples.vcf.gz" I=$inputdir"/Chr6_all_samples.vcf.gz" I=$inputdir"/Chr7_all_samples.vcf.gz" I=$inputdir"/Chr8_all_samples.vcf.gz" \
I=$inputdir"/Chr9_all_samples.vcf.gz" I=$inputdir"/Chr10_all_samples.vcf.gz" I=$inputdir"/Chr11_all_samples.vcf.gz" I=$inputdir"/Chr12_all_samples.vcf.gz" \
I=$inputdir"/Chr13_all_samples.vcf.gz" I=$inputdir"/Chr14_all_samples.vcf.gz" I=$inputdir"/Chr15_all_samples.vcf.gz" I=$inputdir"/Chr16_all_samples.vcf.gz" \
I=$inputdir"/Chr17_all_samples.vcf.gz" I=$inputdir"/Chr18_all_samples.vcf.gz" I=$inputdir"/Chr19_all_samples.vcf.gz" I=$inputdir"/Chr20_all_samples.vcf.gz" \
I=$inputdir"/Chr21_all_samples.vcf.gz" I=$inputdir"/Chr22_all_samples.vcf.gz" I=$inputdir"/Chr23_all_samples.vcf.gz" I=$inputdir"/Chr24_all_samples.vcf.gz" \
O=$inputdir"/"$output


echo "Merging of all chromosomes in a single VCF file is done!\nYou can now run the scripts 4_1_variants_stats.sh for output statistics and 4_2_variants_filtering.sh to filter the VCF file."
