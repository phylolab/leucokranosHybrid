#!/usr/bin/bash -l

#SBATCH --job-name=varstats
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --output=outscripts/%x.%j.out
#SBATCH --error=outscripts/%x.%j.err
#SBATCH --time=3:00:00

#### The script is used to generate statistics from VCF files ####
#### The statistics can be calculated (1) before VCF file filtering, to determine the correct thresholds, and (2) after filtering
#### Author: Wan-Ting Huang


###########################
# Arguments to pass
workindir=${1-/scratch/sschmi13/leucokranos/4_variant_filtering/4_2_filtered_vcf/final_filtered_vcf}			
inputfile=${2-all_chr_no_lowdp_indv_all.gatk_filtered.no_mono.biall.mac.maxmiss.mindp.maxdp.qual.vcf.gz}					

############################################################

####################
# LOAD MODULES
####################

module load vcftools/0.1.16

####################
# BEGINNING OF SCRIPT
####################

cd $workindir

# Get prefix name for output
prefix=${inputfile%.vcf.gz}

# Calculate allele frequency
vcftools --gzvcf $inputfile --freq2 --out $prefix --max-alleles 2
echo "Done calculating allele frequency"

# Calculate mean depth per individual
vcftools --gzvcf $inputfile --depth --out $prefix
echo "Done calculating mean depth for each sample"

# Calculate mean depth per site
vcftools --gzvcf $inputfile --site-mean-depth --out $prefix
echo "Done calculating mean depth for each SNP"

# Calculate site quality
vcftools --gzvcf $inputfile --site-quality --out $prefix
echo "Done calculating site quality"

# Calculate proportion of missing data per individual
vcftools --gzvcf $inputfile --missing-indv --out $prefix
echo "Done calculating the proportion of missing data for each sample"

# Calculate proportion of missing data per site
vcftools --gzvcf $inputfile --missing-site --out $prefix
echo "Done calculating the proportion of missing data for each SNP"

# Calculate heterozygosity and inbreeding coefficient per individual
vcftools --gzvcf $inputfile --het --out $prefix
echo "Done calculating heterozygosity and inbreeding coefficient for each sample"

# Identify location of singletons, and the individual they occur in
vcftools --gzvcf $inputfile --singletons --out $prefix
echo "Done identifying singletons"


####################
# END
####################

				
echo "All global statistics about the VCF files are calculated!"				
	

