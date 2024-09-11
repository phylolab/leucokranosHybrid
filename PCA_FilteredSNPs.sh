#!/usr/bin/bash -l

#SBATCH --job-name=PCA
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --output=outscripts/PCA.out
#SBATCH --error=outscripts/PCA.err
#SBATCH --time=1:00:00

###########################
# Arguments to pass:
workindir=.					# Working directory. Default: Current directory
inputfileDir=4B_FilteredCalling			# Name of the directory with the filtered callset (should be in workindir). Default: 4B_FilteredCalling
outputfileDir=5_PCA_FilteredSNPs			# Name of the directory with files required to do PCA. Default: 5_PCA_FilteredSNPs

############################################################


# Check settings
echo "Using following arguments: 

Working directory:	"$workindir"

Prune:
	-- "$inputfileDir"/GATK_RecStandards/AllSamples_SNP.Filtered.vcf.gz
	-- "$inputfileDir"/SpeciationGenomicsTutorial/AllSamples_SNP.SG.Filtered.vcf.gz

Outputs of pruning and inputs of PCA:
	-- "$outputfileDir"/AllSamples_SNP.Filtered.prune.in
	-- "$outputfileDir"/AllSamples_SNP.SG.Filtered.prune.in

"

echo "
If settings are wrong, please look at instructions and arguments to pass

"

####################
# LOAD MODULES
####################

module load plink/1.9-beta6.27

####################
# BEGINNING OF SCRIPT
####################
cd $workindir

mkdir -p 5_PCA_FilteredSNPs

# Prune the variants that are in linkage
plink --vcf $inputfileDir"/GATK_RecStandards/AllSamples_SNP.Filtered.vcf.gz" \
	--double-id --allow-extra-chr \
	--set-missing-var-ids @:# \
	--indep-pairwise 50 10 0.2 \
	 --out $outputfileDir"/AllSamples_SNP.Filtered"
 

plink --vcf $inputfileDir"/SpeciationGenomicsTutorial/AllSamples_SNP.SG.Filtered.vcf.gz" \
	--double-id --allow-extra-chr \
	--set-missing-var-ids @:# \
	--indep-pairwise 50 10 0.2 \
	 --out $outputfileDir"/AllSamples_SNP.SG.Filtered"


# Perform PCA on pruned datasets
plink --vcf $inputfileDir"/GATK_RecStandards/AllSamples_SNP.Filtered.vcf.gz" \
 --double-id --allow-extra-chr \
 --set-missing-var-ids @:# \
 --extract $outputfileDir"/AllSamples_SNP.Filtered.prune.in" \
 --make-bed --pca \
 --out $outputfileDir"/AllSamples_SNP.Filtered"
 
 
plink --vcf $inputfileDir"/SpeciationGenomicsTutorial/AllSamples_SNP.SG.Filtered.vcf.gz" \
 --double-id --allow-extra-chr \
 --set-missing-var-ids @:# \
 --extract $outputfileDir"/AllSamples_SNP.SG.Filtered.prune.in" \
 --make-bed --pca \
 --out $outputfileDir"/AllSamples_SNP.SG.Filtered"


echo "Script has finished!"









