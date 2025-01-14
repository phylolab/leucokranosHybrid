#!/usr/bin/bash -l

#SBATCH --job-name=rmv_indv
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=outscripts/%x.%j.out
#SBATCH --error=outscripts/%x.%j.err
#SBATCH --time=3:00:00

## optional script to remove individuals from the VCF file to generate different datasets
## the script has to be adapted to your own file structure and sample names 
## author: Sarah Schmid

## load module
module load vcftools/0.1.16

## remove individuals to generate the different datasets

#dataset1 (parental species and leucokranos, without the clarkii outgroup)
vcftools --gzvcf all_chr_all_samples.gatk_filtered.vcf.gz --keep chryso_sand_leuco.list --recode --stdout | gzip -c > all_chr_no_lowdp_indv_chryso_sand_leuco.gatk_filtered.vcf.gz

#dataset2 (all individuals including clarkii as outgroup)
vcftools --gzvcf all_chr_all_samples.gatk_filtered.vcf.gz --remove-indv LU11 --remove-indv LU61 --remove-indv LU24 --recode --stdout | gzip -c > all_chr_no_lowdp_indv_all.gatk_filtered.vcf.gz

## end of script

echo "The two datasets were generated, please proceed to filtering of the variants using the 4_2_variant_filtering.sh script."