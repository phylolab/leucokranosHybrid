#### GATK pipeline steps until Indel Realignment Step
# Authors: Baptiste Micheli, Sarah Schmid


# The input file is a sorted bam file obtained with mapping alghorithm (i.e. BWA)

### usage python GATK_pipeline.py working_directory MAPPING.sorted.bam sample_id max_cpu max_mem reference_genome
### to run it on multiple samples at the same time in Slurm, please use the companion script 3_1_run_gatk_preprocessing.sh

import os
import sys
import subprocess
import time

################################################################################
# FUNCTIONS #####################################################
################################################################################


#Create directory if it doesn't exist yet
def newDir(path,dirname):					#name of the function and variables
	if not os.path.exists(path+dirname):
		os.makedirs(path+dirname)
	return path+dirname+'/'

#Function that check that input file exists
def check_input_file(file_name):
	if not os.path.exists(file_name):
		print("ERROR: input file %s does not exist\n" % (file_name))
		sys.exit()

#Check index for reference
def check_index(file_name):
	if not os.path.exists(file_name):
		return False
	else:
		return True

#Run line (when using a process with another program)
def run_cline(cline):
	"""this first function runs one cline as a subprocess"""
	return subprocess.call(cline, shell=True)

#Get readgroup ID (see https://angus.readthedocs.io/en/2017/GATK_pipeline.html)
def get_RG_parameters(output_file_sam_header_name, sample_ID):
	input_file = open(output_file_sam_header_name, "r")
	for line in input_file:
		line = line.rstrip() 							#rstrip: remove last line of the file (which is useless and empty). L.rstrip("v") : remove first v from the left. 
		hiseq_info = line.split("\t")[0] 				#split header line at each tab, we just want the first position (which is 0 in python)
		hiseq_info_short = hiseq_info.split(":")[0:4]	#split the header at each ":" and keep only the 4 first terms
		RGID = "_".join(hiseq_info_short)				#join with a "_" all terms extracted
	input_file.close()									#important to always close the file
		
	RGLB = sample_ID 									#library ID (e.g. 006_BIC; or we can also add an information about the barcode for example, but not necessary)
	RGPL = "Illumina" 									#platform (e.g. Illumina, Solid)
	RGSM = sample_ID 									#sample ID (e.g. 006_BIC)
	RGPU = RGID + "." + RGLB							#platform unit
	return RGID, RGLB, RGPL, RGSM, RGPU					#exit the terms from the function to use them afterwards



################################################################################
# PIPELINE #####################################################
################################################################################

start=time.time()

work_dir = sys.argv[1]
mapping_file = sys.argv[2]
sample_ID = sys.argv[3]
CPU=sys.argv[4]
max_mem=sys.argv[5]
ref_genome = sys.argv[6]

# /!\ Check if the paths are correct /!\
# We have to copy the picard-tools exe in our working directory, and then change inside the java into the java full path

picard_exe = "/software/UHTS/Analysis/picard-tools/2.9.0/bin/picard-tools"
GATK_exe = "/software/UHTS/Analysis/GenomeAnalysisTK/4.4.0.0/bin/GenomeAnalysisTK"
samtools_exe = "/software/UHTS/Analysis/samtools/0.1.19/bin/samtools"
java_exe = "/software/Development/java/1.8.0_242/bin/java"
bamUtil_exe ="/scratch/sschmi13/leucokranos/2_bwa_mapping/bamUtil/bin/bam"

output_dir = newDir(work_dir, "3_1_gatk_preprocessing")
tmp_dir=newDir(output_dir, "tmp")

################
# Create header file
################

original_bam = mapping_file
output_file_sam_header_name = output_dir + "/" + sample_ID+".header"

check_input_file(original_bam)

if os.path.exists(output_file_sam_header_name):
	print("%s.header already exists.\nThe creation of the header will not be re-performed.\nPlease erase the %s.header file if you wish to reperform this step\n"  % (sample_ID, sample_ID))
else:
	RG_extract_cline = "%s view %s | head -1 > %s" % (samtools_exe, mapping_file, output_file_sam_header_name)
	run_cline(RG_extract_cline)
print("\nHeader created!\n")
	 
################
# CleanSam with picard-tools
################

original_bam = mapping_file
output_file_clean =  output_dir + "/" + sample_ID+".cleaned.bam"
print("################\nCleanSam Step\n################\n")

check_input_file(original_bam)

if os.path.exists(output_file_clean):
	print("%s.cleaned.bam already exists.\nThe step CleanSam will not be re-performed.\nPlease erase the %s.cleaned.bam file if you wish to reperform this step\n"  % (sample_ID, sample_ID))
else:
	CleanSam_cline = "%s CleanSam I=%s O=%s" % (picard_exe, original_bam, output_file_clean)
	run_cline(CleanSam_cline)
print("\nCleanSam finished!\n")


################
# Add read groups
################

input_file_RG =  output_dir + "/" + sample_ID+".cleaned.bam"
output_file_RG = output_dir + "/" + sample_ID+".cleaned.rg.bam"
print("################\nAddOrReplaceReadGroups Step\n################\n")

check_input_file(input_file_RG)				#function created at the beginning to check whether a file is present or not
check_input_file(output_file_sam_header_name)

if os.path.exists(output_file_RG):
	print("%s.cleaned.rg.bam already exists.\nThe step AddOrReplaceReadGroups will not be re-performed.\nPlease erase the %s.cleaned.rg.bam file if you wish to reperform this step\n"  % (sample_ID, sample_ID))
else:
	RGID, RGLB, RGPL, RGSM, RGPU = get_RG_parameters(output_file_sam_header_name,sample_ID)
	RG_cline="%s AddOrReplaceReadGroups I=%s O=%s RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s" % (picard_exe, input_file_RG, output_file_RG, RGID, RGLB, RGPL, RGPU, RGSM)
	run_cline(RG_cline)
print("\nAddOrReplaceReadGroups finished!\n")


################
# Fix paired-end pairs (Fix SAM FLAGs)
################

input_file_FixPairs = output_dir + "/" + sample_ID+".cleaned.rg.bam"
output_file_FixPairs = output_dir + "/" + sample_ID+".cleaned.rg.fixpairs.bam"
print("################\nsamtools fixmate Step\n################\n")

check_input_file(input_file_FixPairs)

if os.path.exists(output_file_FixPairs):
	print("%s.cleaned.rg.fixpairs.bam already exists.\nThe step samtools fixmate will not be re-performed.\nPlease erase the %s.cleaned.rg.fixpairs.bam file if you wish to reperform this step\n"  % (sample_ID, sample_ID))
else:
	fixmate_cline="%s fixmate %s %s" % (samtools_exe, input_file_FixPairs, output_file_FixPairs)
	run_cline(fixmate_cline)
print("\nsamtools fixmate finished!\n")


################
# Remove secondary alignments
################

input_file_NoSecond = output_dir + "/" + sample_ID+".cleaned.rg.fixpairs.bam"
output_file_NoSecond = output_dir + "/" + sample_ID+".cleaned.rg.fixpairs.nosecondary.bam"
print("################\nRemove secondary alignments Step\n################\n")

check_input_file(input_file_NoSecond)

if os.path.exists(output_file_NoSecond):
	print("%s.cleaned.rg.fixpairs.nosecondary.bam already exists.\nThe step 'Remove secondary alignments' will not be re-performed.\nPlease erase the %s.cleaned.rg.fixpairs.nosecondary.bam file if you wish to reperform this step\n" % (sample_ID, sample_ID))
else:
	no_sec_cline = "%s view -bh -F 256 -o %s %s" % (samtools_exe, output_file_NoSecond, input_file_NoSecond)
	run_cline(no_sec_cline)

print("\n'Remove secondary alignments' finished!\n")


################
# Sort the bam file
################
input_file_sort =  output_dir + "/" + sample_ID+".cleaned.rg.fixpairs.nosecondary.bam"
output_file_sort =  output_dir + "/" + sample_ID+".cleaned.rg.fixpairs.nosecondary.sort.bam"
print("################\nSorting Resulting BAM\n################\n")

check_input_file(input_file_sort)

if os.path.exists(output_file_sort):
	print("%s.cleaned.rg.fixpairs.nosecondary.sort.bam already exists.\nThe sorting will not be reperformed.\nPlease erase the %s.cleaned.rg.fixpairs.nosecondary.sort.bam file if you wish to re-sort it\n"  % (sample_ID, sample_ID))
else:
	sorting_cline = "%s SortSam I=%s O=%s SORT_ORDER=coordinate" % (picard_exe, input_file_sort, output_file_sort)
	run_cline(sorting_cline)

print("\nSorting finished!\n")


################
# MARK DUPLICATES
################
input_file_markDup = output_dir + "/" + sample_ID+".cleaned.rg.fixpairs.nosecondary.sort.bam"
output_file_markDup = output_dir + "/" + sample_ID+".dedupped.bam"
output_file_markDup_metric = output_dir + "/" + sample_ID+".dedupped.metrics.txt"
print("################\nMark Duplicated Step\n################\n")

check_input_file(input_file_markDup)

if os.path.exists(output_file_markDup):
	print("%s.dedupped.bam already exists.\nMark Duplicated step will not be reperformed.\nPlease erase the %s.dedupped.bam file if you wish to reperform this step\n"  % (sample_ID, sample_ID))
else:
	makdup_cline = "%s MarkDuplicates I=%s O=%s M=%s  CREATE_INDEX=true" % (picard_exe,input_file_markDup, output_file_markDup, output_file_markDup_metric)
	run_cline(makdup_cline)

print("\nMark Duplicates finished!\n")


################
# Check the presence of Reference for GATK, and if not present generate it
################
check_input_file(ref_genome)
fai_index_file = ref_genome+".fai"
ref_index = ref_genome.rfind(".")
picard_index_file = ref_genome[:ref_index]+".dict"
print("################\nGenerate Reference Index\n################\n")


if check_index(fai_index_file):
	print("Reference.fai file already generating.\nContinuing to following step.\n If you want to regenerate the file, please erase Reference.fai\n")
else:
	fai_index_cline="%s faidx %s" % (samtools_exe, ref_genome)
	run_cline(fai_index_cline)


if check_index(picard_index_file):
	print("Reference.dict file already generating.\nContinuing to following step.\n If you want to regenerate the file, please erase Reference.fai\n")
else:
	pic_index_cline="%s CreateSequenceDictionary R=%s O=%s" % (picard_exe, ref_genome, picard_index_file)
	run_cline(pic_index_cline)

print("\nReference indexing finished!\n")



################
# Reorder BAM file in the same order as the reference
################
input_file_reorder =  output_dir + "/" + sample_ID+".dedupped.bam"
output_file_reorder =  output_dir + "/" + sample_ID+".dedupped.reorder.bam"

print("################\nReorder BAM file\n################\n")

check_input_file(input_file_reorder)

if os.path.exists(output_file_reorder):
	print("%s.dedupped.reorder.bam already exists.\n Reorder BAM file will not be reperformed.\nPlease erase the %s.dedupped.reorder.bam file if you wish to reperform this step\n"  % (sample_ID, sample_ID))
else:
	reorderBAM_cline = "%s ReorderSam INPUT=%s OUTPUT=%s REFERENCE=%s" % (picard_exe, input_file_reorder, output_file_reorder, ref_genome)
	run_cline(reorderBAM_cline)

print("\nReorder BAM file finished!\n")


################
# Generate BAI file
################
input_file_bai =  output_dir + "/" + sample_ID+".dedupped.reorder.bam"
output_file_bai =  output_dir + "/" + sample_ID+".dedupped.reorder.bam.bai"

print("################\nGeneration of BAI file\n################\n")

check_input_file(input_file_bai)

if os.path.exists(output_file_bai):
	print("%s.dedupped.reorder.bam.bai already exists.\n Generation of BAI file will not be reperformed.\nPlease erase the %s.dedupped.reorder.bai file if you wish to reperform this step\n"  % (sample_ID, sample_ID))
else:
	generationBAI_cline = "%s index %s" % (samtools_exe, input_file_bai)
	run_cline(generationBAI_cline)

print("\nGeneration of BAI file finished!\n")


################
# Validation of BAM file
################

input_file_validation =  output_dir + "/" + sample_ID+".dedupped.reorder.bam"
output_file_validation =  output_dir + "/" + sample_ID+".validation.txt"
ref_genome_no_ext = ref_genome.replace(".fna", "")
 
print("################\nValidation of BAM file\n################\n")

check_input_file(input_file_validation)

if os.path.exists(output_file_validation):
	print("%s.validation.txt already exists.\n Validation of BAM file will not be reperformed.\nPlease erase the %s.validation.txt file if you wish to reperform this step\n"  % (sample_ID, sample_ID))
else:
	validateBAM_cline = "%s validate --in %s --refFile %s --so_coord --verbose 2>%s" % (bamUtil_exe, input_file_validation, ref_genome_no_ext, output_file_validation)
	run_cline(validateBAM_cline)

print("\nValidation of BAM file finished!\n")

print("\ngatk_preprocessing script finished!\nPlease run 3_2_index_ref_genome.sh if your reference genome is not indexed for GATK or 3_3_gatk_haplotype_call.sh for haplotype calling\n")

end=time.time()
total_time=end-start
print(total_time)