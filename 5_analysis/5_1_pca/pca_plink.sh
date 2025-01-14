#!/usr/bin/bash -l

#SBATCH --job-name=pca
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=outscripts/%x.%j.out
#SBATCH --error=outscripts/%x.%j.err
#SBATCH --time=1:00:00

#### This script is used to prune the variants and perform a PCA on the pruned dataset ####
#### Author: Wan-Ting Huang

###########################
# Arguments to pass
workindir=${1-.}
inputfile=${2-input_vcf/all_chr_no_lowdp_indv_all.gatk_filtered.no_mono.biall.mac.maxmiss.mindp.maxdp.qual.vcf.gz}
outputdir=${3-pca_output}

############################################################


# Check settings

echo "Using following arguments: 

Working directory:	"$workindir"

Prune:
	-- "$inputfile


####################
# LOAD MODULES
####################

module load plink/1.9-beta6.27 htslib/1.17

####################
# BEGINNING OF SCRIPT
####################

cd $workindir

mkdir -p $outputdir

file_name=${inputfile%.vcf.gz}
prefix=${file_name##*/}

# Prune the variants that are in linkage
plink --vcf $inputfile \
	--double-id --allow-extra-chr \ # additional chromosome beyond human set authorised
	--set-missing-var-ids @:# \     # giving a SNPs ID for each variable (which we set here to chrom:pos)
	--indep-pairwise 50 10 0.2 \    # pruning (50 kb windows, 10 bp window step, 0.2 r2)
	 --out $outputdir"/"$prefix

# Perform PCA on pruned dataset
plink --vcf $inputfile \
 --double-id --allow-extra-chr \
 --set-missing-var-ids @:# \
 --extract $outputdir"/"$prefix".prune.in" \
 --make-bed --pca \
 --out $outputdir"/"$prefix
 
 
echo "Pruning and PCA done!"
