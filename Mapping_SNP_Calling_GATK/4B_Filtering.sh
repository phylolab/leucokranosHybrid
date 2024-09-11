#!/usr/bin/bash -l

#SBATCH --job-name=Filtering
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --output=outscripts/Filter.out
#SBATCH --error=outscripts/Filter.err
#SBATCH --time=2:00:00

#### The script is used to filter the variants ####
# We tried two kinds of filtering methods with some different thresholds. Comment out one of them if you just want to use one method.

###########################
# Arguments to pass:
workindir=.					# Working directory. Default: Current directory
ReferenceGenome=0_AclarkiiReference/AclarkiiGenome.Chr.fna
inputfileDir=3C_CallingOutput			# Name of the directory with the merged gVCF files each sample (should be in workindir). Default: 3C_CallingOutput
outputfileDir=4B_FilteredCalling			# Name of the directory with a database workspace. Default: 4B_FilteredCalling

############################################################


# Check settings
echo "Using following arguments: 

Working directory:	"$workindir"
Input directory (if relative path, from "$workindir") :	"$inputfileDir"
Output directory: "$outputfileDir


echo "
If settings are wrong, please look at instructions and arguments to pass

"

####################
# LOAD MODULES
####################

module load openjdk/17.0.8.1_1 gatk/4.4.0.0 picard/3.0.0 samtools/1.17 bamtools/2.5.2 vcftools/0.1.16

####################
# BEGINNING OF SCRIPT
####################
cd $workindir


######################
# Follow the filtering steps recommended by GATK
######################
# mkdir -p 4B_FilteredCalling/GATK_RecStandards


# # Extract SNPs to an individual file
# echo "Extract only SNPs to a new file"
# 
# gatk SelectVariants \
# 	-R $ReferenceGenome \
# 	-V $inputfileDir"/AllSamples.vcf.gz" \
# 	-select-type SNP \
# 	-O $outputfileDir"/GATK_RecStandards/AllSamples_SNP.vcf.gz"
# 
# echo "-------------"
# 
# # Annotate filter status of the SNPs
# echo "Apply filters to the SNPs with these thresholds:
# -- QD < 2.0
# -- QUAL < 30.0
# -- SOR > 3.0
# -- FS > 60.0
# -- MQ < 40.0
# -- MQRankSum < -12.5
# -- ReadPosRankSum < -8.0
# 
# Please check the thresholds!
# 
# -------------
# 
# Start to annotate the filter status..."
# 
# gatk VariantFiltration \
#     -V $outputfileDir"/GATK_RecStandards/AllSamples_SNP.vcf.gz" \
#     -filter "QD < 2.0" --filter-name "QD2" \
#     -filter "QUAL < 30.0" --filter-name "QUAL30" \
#     -filter "SOR > 3.0" --filter-name "SOR3" \
#     -filter "FS > 60.0" --filter-name "FS60" \
#     -filter "MQ < 40.0" --filter-name "MQ40" \
#     -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#     -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#     -O $outputfileDir"/GATK_RecStandards/AllSamples_SNP.wFilt.vcf.gz"
#     
# echo "Annotation finished!
# -------------
# 
# Now we select those that pass all filters"

# # Subset SNPs that pass all the filters
# gatk SelectVariants \
# 	--exclude-filtered \
# 	-V $outputfileDir"/GATK_RecStandards/AllSamples_SNP.wFilt.vcf.gz" \
# 	-O $outputfileDir"/GATK_RecStandards/AllSamples_SNP.Filtered.vcf.gz"
# 	
# 
# ###################
# # END
# ###################
# 
# echo "Script has finished!"

######################
# Follow the speciation genomics tutorial (SG)
######################
mkdir -p 4B_FilteredCalling/SpeciationGenomicsTutorial

# Set the variables as thresholds
MAF=0.06
NONMISSING=1
QUAL=40
MIN_DEPTH=2
MAX_DEPTH=50

# Filter the variants
echo "Apply filters to the SNPs with these thresholds:
-- MAF = 0.06
-- NONMISSING = 1
-- QUAL = 40
-- MIN_DEPTH = 2
-- MAX_DEPTH = 50

Please check the thresholds!

-------------

Start to filter using VcfTools"

vcftools --gzvcf $inputfileDir"/AllSamples.vcf.gz" \
		--remove-indels \
		--maf $MAF \
		--max-missing $NONMISSING \
		--minQ $QUAL \
		--minDP $MIN_DEPTH \
		--maxDP $MAX_DEPTH \
		--recode --stdout | gzip -c > $outputfileDir"/SpeciationGenomicsTutorial/AllSamples_SNP.SG.Filtered.vcf.gz"


####################
# END
####################


echo "Script has finished!"
