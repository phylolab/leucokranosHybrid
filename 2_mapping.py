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

##################################################
####################  ARG  #######################
##################################################

output_dir = sys.argv[1] #output directory needs to be created beforehand
SampleID = sys.argv[2]
RefGenome = sys.argv[3]
R1 = sys.argv[4]
R2 = sys.argv[5]
CPU = sys.argv[6]

set_R_path =os.system('PATH="/software/R/3.6.1/bin/:$PATH"')
bwa_exe = "/software/UHTS/Aligner/bwa/0.7.17/bin/"
samtools_exe = "/software/UHTS/Analysis/samtools/1.8/bin/"
bamtools_exe = "/software/UHTS/Analysis/bamtools/2.4.1/bin/"
picard_exe = '/software/UHTS/Analysis/picard-tools/2.9.0/bin/'
##################################################
###############  MAPPING PIPELINE  ###############
##################################################

## GENERATE INDEX REFERENCE FOR BWA ##
#print "################\nGenerate Reference BWA Index\n################\n"

#check_input_file(RefGenome)

#amb_index = RefGenome+".amb"
#ann_index = RefGenome+".ann"
#bwt_index = RefGenome+".bwt"
#pac_index = RefGenome+".pac"
#sa_index = RefGenome+".sa"

#if check_index(amb_index) and check_index(ann_index) and check_index(bwt_index) and check_index(pac_index) and check_index(sa_index):
#print "All the index files for BWA are already generating.\nContinuing to following step.\n If you want to regenerate the files, please erase them\n"
#else:
#bwa_index_cline = "%sbwa index %s" % (bwa_exe, RefGenome)
#run_cline(bwa_index_cline)

#print "\nReference indexing finished!\n"

## MAPPING ##
print("################\nBWA Mapping\n################\n")

output_file_BWA =  output_dir + "/" + SampleID+".BWAMapped.sam"

check_input_file(R1)
check_input_file(R2)

if os.path.exists(output_file_BWA):
	print("%s.BWAMapped.sam already exists.\nThe step BWA Mampping will not be re-performed.\nPlease erase the %s.BWAMapped.sam file if you wish to reperform this step\n" % (SampleID, SampleID))
else:
	bwa_map_cline = "%sbwa mem -M -t %s %s %s %s > %s" % (bwa_exe, CPU, RefGenome, R1, R2, output_file_BWA  )
	run_cline(bwa_map_cline)

print("\nBWA Mapping finished!\n")

## Convert SAM TO BAM ##
print("################\nSam to Bam Step\n################\n")

input_file_sam = output_dir + "/" + SampleID+".BWAMapped.sam"
output_file_bam =   output_dir + "/" + SampleID+".BWAMapped.bam"

check_input_file(input_file_sam)

if os.path.exists(output_file_bam):
	print("%s.BWAMapped.bam already exists.\nThe step Sam to Bam will not be re-performed.\nPlease erase the %s.BWAMapped.bam file if you wish to reperform this step\n" % (SampleID, SampleID))
else:
	sam_to_bam_cline = "%ssamtools view -@ %s -bS %s > %s" % (samtools_exe, CPU, input_file_sam, output_file_bam)
	print(sam_to_bam_cline)
	run_cline(sam_to_bam_cline)

print("\nSam to Bam Step finished!\n")

## SORTING BAM FILE ##
print("################\nSorting Bam Step\n################\n")

input_file_bam = output_dir + "/" + SampleID+".BWAMapped.bam"
output_file_sort =   output_dir + "/" + SampleID+".BWAMapped.Sorted.bam"

check_input_file(input_file_bam)

if os.path.exists(output_file_sort):
	print("%s.BWAMapped.Sorted.bam already exists.\nSorting will not be re-performed.\nPlease erase the %s.BWAMapped.Sorted.bam file if you wish to reperform this step\n" % (SampleID, SampleID))
else:
	sorting_cline = "%ssamtools sort -@ %s %s -o %s " % (samtools_exe, CPU, input_file_bam, output_file_sort)
	run_cline(sorting_cline)

print("\nSorting finished!\n")

 ## GENERATE INDEX BAI ##
print("################\nCreating .bai index Step\n################\n")

input_file_bai = output_dir + "/" + SampleID+".BWAMapped.Sorted.bam"
output_file_bai = output_dir + "/" + SampleID+".BWAMapped.Sorted.bam.bai"

check_input_file(input_file_bai)

if os.path.exists(output_file_bai):
    print("%s.BWAMapped.Sorted.bam.bai already exists.\nThe step indexing .bai will not be re-performed.\nPlease erase the %s.BWAMapped.Sorted.bam.bai file if you wish to reperform this step\n" % (SampleID, SampleID))
else:
 	bai_cline = "%ssamtools index %s" % (samtools_exe, input_file_bai)
 	run_cline(bai_cline)

print("\n.bai indexing finished!\n")

## GENERATE STATISTICS ##
print("################\Alignement statistics Step\n################\n")

input_file_bam =  output_dir + "/" + SampleID+".BWAMapped.Sorted.bam"
output_file_stat = output_dir + "/" + SampleID+".BWAMapped.Stats"
output_InsertSizeMetrics = output_dir + "/" + SampleID+".BWAMapped.InsertSizeMetrics.txt"
output_InsertSizeHisto = output_dir + "/" + SampleID+".BWAMapped.InsertSizeHisto.pdf"

check_input_file(input_file_bam)

if os.path.exists(output_file_stat):
    print("%s.BWAMapped.Stats already exists.\nThe step will not be re-performed.\nPlease erase the %s.BWAMapped.Stats file if you wish to reperform this step\n" % (SampleID, SampleID))
else:
    stat_cline = "%sbamtools stats -in %s > %s " % (bamtools_exe, input_file_bam, output_file_stat)
    run_cline(stat_cline)

if os.path.exists(output_InsertSizeMetrics) and os.path.exists(output_InsertSizeMetrics):
    print("%s.BWAMapped.InsertSizeMetrics.txt and %s.BWAMapped.InsertSizeHisto.pdf and already exist.\nThe step will not be re-performed.\nPlease erase %s.BWAMapped.InsertSizeMetrics.txt and %s.BWAMapped.InsertSizeHisto.pdf files if you wish to reperform this step\n" % (SampleID, SampleID,SampleID,SampleID))
else:
    stat_insertsize = "%spicard-tools CollectInsertSizeMetrics I=%s O=%s H=%s" % (picard_exe, input_file_bam, output_InsertSizeMetrics, output_InsertSizeHisto)
    run_cline(stat_insertsize)

print("################\nAlignment Statistics computation finishes!\n################\n")
print("\nMapping script finished!\nNext step, please run Processing.py\n")

end=time.time()
total_time=end-start
print(total_time)
