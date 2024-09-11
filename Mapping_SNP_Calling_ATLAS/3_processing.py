## Pipeline to generate vcf file from a BAM file (python 3 adapted)

## Overview

## 1) Trimming reads and remove adapters 1_trimming.py
## 2) Mapping trimmed reads to references Fasta file (BWA) 2_mapping.py
## 3) Processing mapping file + STATS (Samtools - Picard -ATLAS -bamtools) 3_processing.py
## 4) Variant Call (ATLAS) 4_glf.py

## Software version (september 2021)

##  - BWA : 0.7.17
##  - Samtools : 1.8
##  - Picard-tools : 2.24.0
##  - ATLAS : 0.99
##  - bamtools : 2.4.1

##################################################
####################  IMPORT  ####################
##################################################

import os
import sys
import subprocess
import time

start=time.time()
##################################################
####################  FUNCTIONS  #################
##################################################

## create directory if it doesn't exist yet
def newDir(path,dirname):
	if not os.path.exists(path+dirname):
		os.makedirs(path+dirname)
	return path+dirname+'/'

## function that check that input file exists
def check_input_file(file_name):
	if not os.path.exists(file_name):
		print("ERROR: input file %s does not exist\n" % (file_name))
		sys.exit()

## run line
def run_cline(cline):
	"""this first function runs one cline as a subprocess"""
	return subprocess.call(cline, shell=True)

#Check index for reference
def check_index(file_name):
	if not os.path.exists(file_name):
		return False
	else:
		return True

#Get readgroup ID
def get_RG_parameters(output_file_sam_header_name, SampleID):
    input_file = open(output_file_sam_header_name, "r")
    for line in input_file:
        line = line.rstrip()                             #rstrip: remove last line of the file (which is useless and empty). L.rstrip("v") : remove first v from the left.
        hiseq_info = line.split("\t")[0]                 #split header line at each tab, we just want the first position (which is 0 in python)
        hiseq_info_short = hiseq_info.split(":")[0:4]    #split the header at each ":" and keep only the 4 first terms
        RGID = "_".join(hiseq_info_short)                #join with a "_" all terms extracted
    input_file.close()                                    #important to always close the file

    RGLB = SampleID                                     #library ID (e.g. 006_BIC; or we can also add an information about the barcode for example, but not necessary)
    RGPL = "Illumina"                                     #platform (e.g. Illumina, Solid)
    RGSM = SampleID                                     #sample ID (e.g. 006_BIC)
    RGPU = RGID + "." + RGLB                            #platform unit
    return RGID, RGLB, RGPL, RGSM, RGPU                    #exit the terms from the function to use them afterwards

##################################################
####################  ARG  #######################
##################################################

WD = sys.argv[1] #important to have the `/` at the end of the output path
SampleID = sys.argv[2]
mapping_file = sys.argv[3]
CPU = sys.argv[4]

output_dir = newDir(WD, "processing_out")
set_R_path =os.system('PATH="/software/R/3.6.1/bin/:$PATH"')
samtools_exe = "/software/UHTS/Analysis/samtools/1.8/bin/samtools"
bamtools_exe = "/software/UHTS/Analysis/bamtools/2.4.1/bin/"
picard_exe = 'picard'
atlas_exe = "singularity run /dcsrsoft/singularity/containers/atlas-0.99.sif"
# when working on Curnagl, add the following two lines in the sbatch script to run atlas
# export SINGULARITY_BINDPATH="/users,/work,/scratch"
# module load singularity


##################################################
#############  Processing PIPELINE  ##############
##################################################


## Generate Header ##
print("################\nCreate Header Step\n################\n")
original_bam = mapping_file
output_file_sam_header_name = output_dir + "/" + SampleID+".BWAMapped.header"

check_input_file(original_bam)

if os.path.exists(output_file_sam_header_name):
    print("%s.BWAMapped.header already exists.\nThe creation of the header will not be re-performed.\nPlease erase the %s.BWAMapped.header file if you wish to reperform this step\n"  % (SampleID, SampleID))
else:
    RG_extract_cline = "%s view %s | head -1 > %s" % (samtools_exe, mapping_file, output_file_sam_header_name)
    run_cline(RG_extract_cline)

print("\nHeader created!\n")

## Add Read Groups ##
print("################\nAddOrReplaceReadGroups Step\n################\n")
input_file_RG =  mapping_file
output_file_RG = output_dir + "/" + SampleID+".BWAMapped.Sorted.RG.bam"

check_input_file(input_file_RG)
check_input_file(output_file_sam_header_name)

if os.path.exists(output_file_RG):
    print("%s.BWAMapped.Sorted.RG.bam already exists.\nThe step AddOrReplaceReadGroups will not be re-performed.\nPlease erase the %s.BWAMapped.Sorted.RG.bam file if you wish to reperform this step\n"  % (SampleID, SampleID))
else:
    RGID, RGLB, RGPL, RGSM, RGPU = get_RG_parameters(output_file_sam_header_name,SampleID)
    RG_cline="%s AddOrReplaceReadGroups I=%s O=%s RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s" % (picard_exe, input_file_RG, output_file_RG, RGID, RGLB, RGPL, RGPU, RGSM)
    run_cline(RG_cline)

print("\nAddOrReplaceReadGroups finished!\n")

 ## Generate Index Bai ##
print("################\nCreating .bai index Step\n################\n")

input_file_bai = output_dir + "/" + SampleID+".BWAMapped.Sorted.RG.bam"
output_file_bai = output_dir + "/" + SampleID+".BWAMapped.Sorted.RG.bam.bai"

check_input_file(input_file_bai)

if os.path.exists(output_file_bai):
    print("%s.BWAMapped.Sorted.RG.bam.bai already exists.\nThe step indexing .bai will not be re-performed.\nPlease erase the %s.BWAMapped.Sorted.RG.bam.bai file if you wish to reperform this step\n" % (SampleID, SampleID))
else:
 	bai_cline = "%s index %s" % (samtools_exe, input_file_bai)
 	run_cline(bai_cline)

print("\nindex finished!\n")

## Filter BAM ##
print("################\nFiltering Step\n################\n")
input_file_filter = output_dir + "/" + SampleID+".BWAMapped.Sorted.RG.bam"
output_file_filter = output_dir + "/" + SampleID+".BWAMapped.filtered.bam"

check_input_file(input_file_filter)

if os.path.exists(output_file_filter):
	print("%s.BWAMapped.filtered.bam already exists.\nThe step indexing .bai will not be re-performed.\nPlease erase the %s.BWAMapped.filtered.bam file if you wish to reperform this step\n" % (SampleID, SampleID))
else:
	filter_cline = "%s view -b -f 2 -F 256 -q 30 -@ %s %s > %s" % (samtools_exe, CPU, input_file_filter,output_file_filter)
	run_cline(filter_cline)

print("\nFiltering finished!\n")

 ## Generate Index Bai ##
print("################\nCreating .bai index Step\n################\n")
input_file_bai2 = output_dir + "/" + SampleID+".BWAMapped.filtered.bam"
output_file_bai2 = output_dir + "/" + SampleID+".BWAMapped.filtered.bam.bai"

check_input_file(input_file_bai2)

if os.path.exists(output_file_bai2):
    print("%s.BWAMapped.filtered.bam.bai already exists.\nThe step indexing .bai will not be re-performed.\nPlease erase the %s.BWAMapped.filtered.bam.bai file if you wish to reperform this step\n" % (SampleID, SampleID))
else:
 	bai_cline2 = "%s index %s" % (samtools_exe, input_file_bai2)
 	run_cline(bai_cline2)

print("\nindex finished!\n")

## readOverlap ##
print("################\nCheck Reads Overlap step\n################\n")
input_file_Overlap = output_dir + "/" + SampleID+".BWAMapped.filtered.bam"
output_file_Overlap =  output_dir + "/" + SampleID+".BWAMapped.filtered_overlapStats.txt"
out_file_Overlap = output_dir + "/" + SampleID+".BWAMapped.filtered"

check_input_file(input_file_Overlap)

if os.path.exists(output_file_Overlap):
	print("%s.BWAMapped.filtered_overlapStats.txt already exists.\nThe step Check Reads Overlap will not be re-performed.\nPlease erase the %s.BWAMapped.filtered_overlapStats.txt file if you wish to reperform this step\n" % (SampleID, SampleID))
else:
	Atlas_Overlap_cline = "%s task=readOverlap bam=%s out=%s" % (atlas_exe, input_file_Overlap,out_file_Overlap)
	run_cline(Atlas_Overlap_cline)

print("\nCheck Reads Overlap finished!\n")

## Remove Overlap Reads ##
print("################\nMerge Read Groups Overlap step\n################\n")
input_file_Merge = output_dir + "/" + SampleID+".BWAMapped.filtered.bam"
output_file_Merge =  output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads.bam"
output_file_Merge_Ignore =  output_dir + "/" + SampleID+ ".BWAMapped.filtered_ignoredReads.txt.gz"
out_file_Merge = output_dir + "/" + SampleID + ".BWAMapped.filtered"

check_input_file(input_file_Merge)

if os.path.exists(output_file_Merge) and os.path.exists(output_file_Merge_Ignore):
	print("%s.BWAMapped.filtered_mergedReads.bam and %s.BWAMapped.filtered_ignoredReads.txt.gz already exists.\nThe step Check Reads Overlap will not be re-performed.\nPlease erase the %s.BWAMapped.filtered_mergedRGs.bam and %s.BWAMapped.filtered_ignoredReads.txt.gz file if you wish to reperform this step\n" % (SampleID, SampleID,SampleID,SampleID))
else:
	Atlas_Merge_cline = "%s task=mergeReads bam=%s out=%s" % (atlas_exe, input_file_Merge,out_file_Merge)
	run_cline(Atlas_Merge_cline)

print("\nMerge Reads finished!\n")

## readOverlap AFTER##
print("################\nCheck Reads Overlap step 2 \n################\n")
input_file_Overlap2 = output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads.bam"
output_file_Overlap2 =  output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads_overlapStats.txt"
out_file_Overlap2 = output_dir + "/" + SampleID + ".BWAMapped.filtered_mergedReads"

check_input_file(input_file_Overlap2)

if os.path.exists(output_file_Overlap2):
	print("%s.BWAMapped.filtered_mergedRGs_overlapStats.txt already exists.\nThe step Check Reads Overlap will not be re-performed.\nPlease erase the %s.BWAMapped.filtered_mergedRGs_overlapStats.txt file if you wish to reperform this step\n" % (SampleID, SampleID))
else:
	Atlas_Overlap_cline2 = "%s task=readOverlap bam=%s out=%s" % (atlas_exe, input_file_Overlap2,out_file_Overlap2)
	run_cline(Atlas_Overlap_cline2)

print("\nCheck Reads Overlap after finished!\n")

## Create Depth Mask ##
print("################\nDepth Mask step\n################\n")
input_file_DepthM = output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads.bam"
output_file_DepthM =  output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads_minDepth0_maxDepth40_depthMask.bed"
out_file_DepthM = output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads"

check_input_file(input_file_DepthM)

if os.path.exists(output_file_DepthM):
	print("%s.BWAMapped.filtered_mergedRGs_minDepth0_maxDepth40_depthMask.bed already exists.\nThe step Check Reads Overlap will not be re-performed.\nPlease erase the %s.BWAMapped.filtered_mergedRGs_minDepth0_maxDepth10000000_depthMask.bed file if you wish to reperform this step\n" % (SampleID, SampleID))
else:
	Atlas_DepthM_cline = "%s task=createDepthMask minDepthForMask=0 maxDepthForMask=40 bam=%s out=%s" % (atlas_exe, input_file_DepthM,out_file_DepthM)
	run_cline(Atlas_DepthM_cline)

print("\nDepth Mask finished!\n")

## assess SoftClipping ##
print("################\nAsses SoftClipping step\n################\n")
input_file_SoftCl = output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads.bam"
output_file_SoftCl =  output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads_clippingStats.txt.gz"
out_file_SoftCl = output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads"

check_input_file(input_file_SoftCl)

if os.path.exists(output_file_SoftCl):
	print("%s.BWAMapped.filtered_mergedRGs_clippingStats.txt.g already exists.\nThe step Check Reads Overlap will not be re-performed.\nPlease erase the %s.BWAMapped.filtered_mergedRGs_clippingStats.txt.g file if you wish to reperform this step\n" % (SampleID, SampleID))
else:
	Atlas_SoftCl_cline = "%s task=assessSoftClipping bam=%s out=%s" % (atlas_exe, input_file_SoftCl,out_file_SoftCl)
	run_cline(Atlas_SoftCl_cline)

print("\nAsses SoftClipping finished!\n")

## GENERATE STATISTICS ##
print("################\Alignement statistics Step\n################\n")

input_file_bam =  output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads.bam"
output_file_stat = output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads.Stats"
output_InsertSizeMetrics = output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads.InsertSizeMetrics.txt"
output_InsertSizeHisto = output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads.InsertSizeHisto.pdf"

check_input_file(input_file_bam)

if os.path.exists(output_file_stat):
    print("%s.BWAMapped.filtered_mergedReads.Stats already exists.\nThe step will not be re-performed.\nPlease erase the %s.BWAMapped.filtered_mergedReads.Stats file if you wish to reperform this step\n" % (SampleID, SampleID))
else:
    stat_cline = "%sbamtools stats -in %s > %s " % (bamtools_exe, input_file_bam, output_file_stat)
    run_cline(stat_cline)

if os.path.exists(output_InsertSizeMetrics) and os.path.exists(output_InsertSizeMetrics):
    print("%s.BWAMapped.filtered_mergedReads.InsertSizeMetrics.txt and %s.BWAMapped.filtered_mergedReads.InsertSizeHisto.pdf and already exist.\nThe step will not be re-performed.\nPlease erase %s.BWAMapped.filtered_mergedReads.InsertSizeMetrics.txt and %s.BWAMapped.filtered_mergedReads.InsertSizeHisto.pdf files if you wish to reperform this step\n" % (SampleID, SampleID,SampleID,SampleID))
else:
    stat_insertsize = "%s CollectInsertSizeMetrics I=%s O=%s H=%s" % (picard_exe, input_file_bam, output_InsertSizeMetrics, output_InsertSizeHisto)
    run_cline(stat_insertsize)

print("\nBAM Stats finished!\n")


## GENERATE STATISTICS ATLAS ##
print("################\nStats ATLAS Step\n################\n")

input_file_StatsAtlas = output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads.bam"
output_file_StatsAtlas =  output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads_approximateDepth.txt"
output_file_StatsAtlas2 =  output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads_MQ.txt"
output_file_StatsAtlas3 =  output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads_readLength.txt"
output_file_StatsAtlas4 =  output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads_fragmentStats.txt"
out_file_StatsAtlas = output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads"

check_input_file(input_file_StatsAtlas)

if os.path.exists(output_file_StatsAtlas) and os.path.exists(output_file_StatsAtlas2) and os.path.exists(output_file_StatsAtlas3) and os.path.exists(output_file_StatsAtlas4):
    print("%s.BWAMapped.filtered_mergedReads_ATLAS_STATS.files already exists.\nThe step will not be re-performed.\nPlease erase the %s.BWAMapped.filtered_mergedReads_ATLAS_STATS.files files if you wish to reperform this step\n" % (SampleID, SampleID))
else:
	Atlas_Stats_cline = "%s task=BAMDiagnostics bam=%s out=%s" % (atlas_exe, input_file_StatsAtlas,out_file_StatsAtlas)
	run_cline(Atlas_Stats_cline)


print("\nATLAS Stats finished!\n")

## Pileup Step ##
print("################\nPileup step\n################\n")
input_file_pileup = output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads.bam"
output_file_pileup =  output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads_pileup.txt"
out_file_pileup = output_dir + "/" + SampleID+".BWAMapped.filtered_mergedReads"

check_input_file(input_file_pileup)

if os.path.exists(output_file_pileup):
	print("%s.BWAMapped.filtered_mergedRGs_pileup.txt already exists.\nThe step Check Reads Overlap will not be re-performed.\nPlease erase the %s.BWAMapped.filtered_mergedRGs_pileup.txt file if you wish to reperform this step\n" % (SampleID, SampleID))
else:
	Atlas_pileup_cline = "%s task=pileup bam=%s out=%s" % (atlas_exe, input_file_pileup,out_file_pileup)
	run_cline(Atlas_pileup_cline)

print("\npileup Step finished!\n")

## Organise outputs ##
newDir(output_dir,"Stats_ReadOverlap")
newDir(output_dir,"Stats_DepthMask")
newDir(output_dir,"Stats_SoftClipping")
newDir(output_dir,"Stats_BAMDiagnostics")
newDir(output_dir,"Stats_Pileup")
newDir(output_dir,"Stats_bamtools")

mv_cline="mv processing_out/%s.BWAMapped.filtered_mergedReads_minDepth0_maxDepth40_depthMask.bed processing_out/Stats_DepthMask/" % (SampleID)
run_cline(mv_cline)
mv_cline2="mv processing_out/%s.BWAMapped.filtered_mergedReads_overlapStats.txt processing_out/Stats_ReadOverlap/ ; mv processing_out/%s.BWAMapped.filtered_mergedReads_ignoredReads.txt.gz processing_out/Stats_ReadOverlap" % (SampleID,SampleID)
run_cline(mv_cline2)
mv_cline3="mv processing_out/%s.BWAMapped.filtered_mergedReads_clippingStats.txt.gz processing_out/Stats_SoftClipping/" % (SampleID)
run_cline(mv_cline3)
mv_cline4="mv processing_out/%s.BWAMapped.filtered_mergedReads_approximateDepth.txt processing_out/Stats_BAMDiagnostics ; mv processing_out/%s.BWAMapped.filtered_mergedReads_MQ.txt processing_out/Stats_BAMDiagnostics ; mv processing_out/%s.BWAMapped.filtered_mergedReads_readLength.txt processing_out/Stats_BAMDiagnostics ; mv processing_out/%s.BWAMapped.filtered_mergedReads_fragmentStats.txt processing_out/Stats_BAMDiagnostics " % (SampleID,SampleID,SampleID,SampleID)
run_cline(mv_cline4)
mv_cline5="mv processing_out/%s.BWAMapped.filtered_mergedReads_pileup.txt processing_out/Stats_Pileup/" % (SampleID)
run_cline(mv_cline5)
mv_cline6="mv processing_out/%s.BWAMapped.filtered_mergedReads.Stats processing_out/Stats_bamtools/; mv processing_out/%s.filtered_mergedReads.InsertSizeMetrics.txt processing_out/Stats_bamtools/ ; mv processing_out/%s.filtered_mergedReads.InsertSizeHisto.pdf processing_out/Stats_bamtools/" % (SampleID,SampleID,SampleID)
run_cline(mv_cline6)

print("\Processing script finished!\nNext step, please run SNP_ATLAS.py\n")

end=time.time()
total_time=end-start
print(total_time)
