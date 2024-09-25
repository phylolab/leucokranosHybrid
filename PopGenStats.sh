#!/usr/bin/bash -l

#SBATCH --job-name=PopGenStats
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --output=outscripts/PopGenStats.out
#SBATCH --error=outscripts/PopGenStats.err
#SBATCH --time=48:00:00

#### The script is used to run the scripts of Simon Martin (https://github.com/simonhmartin/genomics_general) ####

##############################
# Argument to pass: 
workdir=.
pyconvertvcf=5_PopGenStats/genomics_general-master/VCF_processing/parseVCF.py
pywindow=5_PopGenStats/genomics_general-master/popgenWindows.py
pyabba=5_PopGenStats/genomics_general-master/ABBABABAwindows.py
vcffile=4B_FilteredCalling/SpeciationGenomicsTutorial/AllSamples_SNP.SG.Filtered.vcf.gz
inputfileDir=5_PopGenStats
inputfile=AllSamples_SNP.SG.Filtered.geno.gz
outputfileDir=5_PopGenStats
poplabel=5_PopGenStats/PopLabel.txt
p1="SAN_AUS"
p2="SAN_PNG"
p3="CHR_PNG"
outgroup="CLA"
minSNPs=100
windowsize=50000

#######################################
# Check the settings

echo "Using following setteings:

Input file:		"$inputfile"

Window type:	Coordinate
Examined window size:		"$windowsize" bp
Windows contain a minimum of "$minSNPs" SNPs

Please check the settings carefully!

-------------------------------------"



####################
# BEGINNING OF SCRIPT
####################
cd $workdir


# Convert vcf.gz file to geno.gz file which is the required format
python $pyconvertvcf -i $vcffile | gzip > $inputfileDir"/"$inputfile


# Compute standard population genomic statistics: pi, Fst, and dxy, for each pair of species
# Each population is listed using -p option
echo "Computing pi, Fst and dxy..."

python $pywindow -f phased --windType coordinate \
		-w $windowsize -m $minSNPs -T 2 \
		-g $inputfileDir"/"$inputfile \
		-o $outputfileDir"/AllSpecies_popgenWindows.w50m100.csv.gz" \
		-p CHR_SLM -p CHR_PNG -p CHR_FIJI -p CHR_PLU \
		-p leuco_BC	-p leuco_F1 -p leuco_F2 -p leuco_BC_chryso \
		-p SAN_PNG_admix -p SAN_PNG -p SAN_SLM -p SAN_AUS -p SAN_miss \
		-p CLA \
		--popsFile $poplabel

echo "Output file:	"$outputfileDir"/AllSpecies_popgenWindows.w50m100.csv.gz

-------------------------------------

Compute Fd between the populations:
	-P1:		A. sandaracinos from Australia (labeled as "$p1")
	-P2:		A. sandaracinos from Papua New Guinea (labeled as "$p2")
	-P3:		A. chrysopterus from Papua New Guinea (labeled as "$p3")
	
	Outgroup:		A. clarkii (labeled as "$outgroup")
	
	"


# Compute Fd between certain populations we're interested in
python $pyabba -f phased \
	-w $windowsize -m $minSNPs -T 2 \
	-g $inputfileDir"/"$inputfile \
	-o $outputfileDir"/ABBABABA.SAN_AUS.SAN_PNG.CHR_PNG.CLA.w50m100.csv.gz" \
	-P1 $p1 -P2 $p2 -P3 $p3 -O $outgroup \
	--popsFile $poplabel
	
	
echo "Script has finished!"


