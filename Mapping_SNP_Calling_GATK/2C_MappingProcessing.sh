#!/usr/bin/bash -l

#SBATCH --job-name=Processing
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --output=outscripts/Processing_%a.out
#SBATCH --error=outscripts/Processing_%a.err
#SBATCH --time=20:00:00

#### The script is used to filter the mapped reads and generate statistics of mapped reads ####


##############################
# Argument to pass: 
SampleList=${1-AllSamples.txt}			# Name of the file with sample names. Default: AllSamples.txt
workindir=${2-.}					# Working directory. Default: Current Directory
inputfileDir=${3-2B_MappingOutput}			# Name of the directory with the input mapping files (should be in workindir). Default: 2B_MappingOutput
QualityThreshold=${4-30}				# Mapping quality threshold. Default 30
AtlasExe=${5-atlas}	#Atlas executable. Default 
##############################


############################################################
## Relaunch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $SampleList | awk '{print $1}')


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $SampleList $workindir $inputfileDir $QualityThreshold $AtlasExe
fi

############################################################
# Check settings
echo "Using following arguments: 

List of Samples:	"$SampleList"
Working directory:	"$workindir"
Mapping Files directory (if relative path, from "$workindir") :	"$inputfileDir"
Mapping quality threshold:	"$QualityThreshold"
Atlas executable: 	"$AtlasExe

echo "
If settings are wrong, please look at instructions and arguments to pass

"
####################
# LOAD MODULES
####################

module load bwa/0.7.17 python/3.11.6 samtools/1.17 picard/3.0.0 bamtools/2.5.2 r/4.3.2

module load miniconda3/22.11.1
conda activate /work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/Conda/envs/ATLAS

####################
# BEGINNNING OF SCRIPT
####################

# Get name of samples for array
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SampleList)

# Print the name of the sample
echo "Working on sample:	"$sample

# Change to working directory
cd $workindir

####################
# Generate file with read name to later extract read group
####################

mkdir -p temp

samtools view "./"$inputfileDir"/"$sample".BWA.Aclarkii.Sorted.bam" | head -n 1 > "./temp/"$sample".header.txt"

####################
# Get information for Readgroups
####################

string="$(cut -f1 "./temp/"$sample".header.txt")"
RGID="$(cut -d':' -f1-4 <<< $string)"
RGLB=$sample
RGPL="Illumina"
RGPU=$RGID"."$RGLB
RGSM=$sample

 
####################
# Add Read Groups
####################
mkdir -p 2C_FilteredMapping


picard AddOrReplaceReadGroups -I "./"$inputfileDir"/"$sample".BWA.Aclarkii.Sorted.bam" \
-O "./2C_FilteredMapping/"$sample".BWA.Aclarkii.Sorted.RG.bam" \
-RGID "$RGID" -RGLB "$RGLB" -RGPL "$RGPL" -RGPU "$RGPU" -RGSM "$RGSM"


####################
# Index the bam file with read group
####################

samtools index "./2C_FilteredMapping/"$sample".BWA.Aclarkii.Sorted.RG.bam"

####################
# Filter BAM 
# Keep only files with FLAGS 2 (read mapped in proper pair) and remove files with  FLAGS 256 (not a primary alignment)
####################

samtools view -b -f 2 -F 256 -q $QualityThreshold \
-@ 2 "./2C_FilteredMapping/"$sample".BWA.Aclarkii.Sorted.RG.bam" \
> "./2C_FilteredMapping/"$sample".BWA.Aclarkii.Sort.Filt.bam"

####################
# Index the filtered bam
####################

samtools index "./2C_FilteredMapping/"$sample".BWA.Aclarkii.Sort.Filt.bam"

####################
# Run Atlas readOverlap
# Generate statistics on overlapping read pairs with ATLAS.
####################

mkdir -p 2C_ProcessingStatistics

$AtlasExe task=readOverlap bam="./2C_FilteredMapping/"$sample".BWA.Aclarkii.Sort.Filt.bam" \
out="./2C_ProcessingStatistics/"$sample".BWA.Aclarkii.Sort.Filt"

####################
# Run atlas_readMerge
#  Generate a new mapping file by setting the quality of one of the overlapping segment paired reads to zero with atlas.
####################

$AtlasExe task=mergeReads bam="./2C_FilteredMapping/"$sample".BWA.Aclarkii.Sort.Filt.bam" \
out="./2C_FilteredMapping/"$sample".BWA.Aclarkii.Sort.Filt"


####################
# Generate Statistics with bamtools
####################

bamtools stats -in "./2C_FilteredMapping/"$sample".BWA.Aclarkii.Sort.Filt_mergedReads.bam" \
> "./2C_ProcessingStatistics/"$sample".BWA.Aclarkii.Sort.Filt_mergedReads.Stats"


####################
# Generate additional mapping statistics (insert size)
####################

picard CollectInsertSizeMetrics -I "./2C_FilteredMapping/"$sample".BWA.Aclarkii.Sort.Filt_mergedReads.bam" \
-O "./2C_ProcessingStatistics/"$sample".Filt.MergedReads.InsertMetrics.txt" \
-H "./2C_ProcessingStatistics/"$sample".Filt.MergedReads.InsertPlot.txt"


####################
# Generate additional mapping statistics wth ATLAS (insert size)
####################

$AtlasExe task=BAMDiagnostics \
bam="./2C_FilteredMapping/"$sample".BWA.Aclarkii.Sort.Filt_mergedReads.bam" \
out="./2C_ProcessingStatistics/"$sample".BWA.Aclarkii.Sort.Filt_mergedReads"


####################
# END
####################

echo "The script for sample "$sample" has finished"
echo "If you wanna run ATLAS createDepthMask, assessSoftClipping or pileup, please look at the end of the 3_MappingProcessing.sh script"

#######################################################################################
#######################################################################################

# following task are not performed by default, as they are computationally intensive and the output is not essential, as they are mostly for statistics purposes. A filter on depth will be done on the VCF in later steps. Softclipping is not considered in callers. If you still want to run the following task, please run it separately

#######################################################################################
#######################################################################################


####################
# Run atlas_depth_mask
# Create a depth mask in BED format
# If still want to run, please check and set the parameters minDepthForMasking and maxDepthForMask

#minDepthForMasking=0
#maxDepthForMask=40
#$AtlasExe task=createDepthMask minDepthForMask=$minDepthForMasking maxDepthForMask=$maxDepthForMasking bam="./FilteredMapping/"$sample".BWA.Aclarkii.Sort.Filt_mergedReads.bam" out="./FilteredMapping/"$sample".BWA.Aclarkii.Sort.Filt_mergedReads"


####################
# Run Assess soft clipping
# Assess soft clipping bases
# Soft-clipped bases are not considered by callers


#$AtlasExe task=assessSoftClipping bam="./FilteredMapping/"$sample".BWA.Aclarkii.Sort.Filt_mergedReads.bam" out="./ProcessingStatistics/"$sample".BWA.Aclarkii.Sort.Filt_mergedReads"



####################
# Generate additional mapping statistics wth ATLAS (pileup)
# !!!!! By default, we do not run it as much memory and computation !!!!!


#$AtlasExe task=pileup bam="./FilteredMapping/"$sample"BWA.Aclarkii.Sort.Filt_mergedReads.bam" out="./ProcessingStatistics/"$sample"BWA.Aclarkii.Sort.Filt_mergedReads"










