#!/usr/bin/bash -l

#SBATCH --job-name=pop_gen_stats
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --output=outscripts/%x.out
#SBATCH --error=outscripts/%x.err
#SBATCH --time=48:00:00

#### The script is used to run the scripts of Simon Martin (https://github.com/simonhmartin/genomics_general) ####
#### It calculate standard population genomics statistics in sliding windows as well as Fd values for introgression detection ####
#### Authors: Wan-Ting Huang, Sarah Schmid

##############################
# Argument to pass: 
workdir=${1-/scratch/sschmi13/6_pop_gen_fd}
pyconvertvcf=${2-genomics_general/VCF_processing/parseVCF.py}
pywindow=${3-genomics_general/popgenWindows.py}
pyabba=${4-genomics_general/ABBABABAwindows.py}
vcffile=${5-../final_vcf/all_chr_no_lowdp_indv_all.gatk_filtered.no_mono.biall.mac.maxmiss.mindp.maxdp.qual.annotate.vcf} #should be UTF-8 encoded, otherwise it won't work
inputdir=${6-genofile}
inputfile=${7-all_chr_no_lowdp_indv_all.gatk_filtered.no_mono.biall.mac.maxmiss.mindp.maxdp.qual.geno.gz}
outputdir=${8-results}
poplabel=${9-pop_label.txt}
p1=${10-"SAN_AUS"}
p2=${11-"SAN_PNG"}
p3=${12-"CHR_PNG"}
p1_2=${10-"CHR_FIJI"}
p2_2=${11-"CHR_PNG"}
p3_2=${12-"SAN_PNG"}
outgroup=${13-"CLA"}
minSNPs=${14-100}
windowsize=${15-50000}

#######################################
# Check the settings

echo "Using following settings:

Input vcf file:		"$vcffile"
Input geno file:		"$inputfile"

Window type:	Coordinate
Examined window size:		"$windowsize" bp
Windows contain a minimum of "$minSNPs" SNPs

Please check the settings carefully!

-------------------------------------"



####################
# BEGINNING OF SCRIPT
####################
cd $workdir

mkdir -p $outputdir
mkdir -p $inputdir



# Convert vcf.gz file to geno.gz file which is the required format
python $pyconvertvcf -i $vcffile | gzip > $inputdir"/"$inputfile


# Compute standard population genomic statistics: pi, Fst, and dxy, for each pair of species
# Each population is listed using -p option
echo "Computing pi, Fst and dxy..."

python $pywindow -f phased --windType coordinate \
		-w $windowsize -m $minSNPs -T 2 \
		-g $inputdir"/"$inputfile \
		-o $outputdir"/all_species_popgen.w50m100.csv.gz" \
		-p CHR_SLM -p CHR_PNG -p CHR_FIJI -p CHR_PLU \
		-p leuco_BC	-p leuco_F1 -p leuco_F2 -p leuco_BC_chryso \
		-p SAN_PNG_admix -p SAN_PNG -p SAN_SLM -p SAN_AUS -p SAN_miss \
		-p CLA \
		--popsFile $poplabel

echo "Output file:	"$outputdir"/all_species_popgen.w50m100.csv.gz

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
	-g $inputdir"/"$inputfile \
	-o $outputdir"/fd.san_aus.san_png.chr_png.cla.w50m100.csv.gz" \
	-P1 $p1 -P2 $p2 -P3 $p3 -O $outgroup \
	--popsFile $poplabel
	
echo "Output file:	"$outputdir"/fd.san_aus.san_png.chr_png.cla.w50m100.csv.gz

-------------------------------------

Compute Fd between the populations:
	-P1:		A. chrysopterus from Fiji (labeled as "$p1_2")
	-P2:		A. chrysopterus from Papua New Guinea (labeled as "$p2_2")
	-P3:		A. sandaracinos from Papua New Guinea (labeled as "$p3_2")
	
	Outgroup:		A. clarkii (labeled as "$outgroup")
	
	"


# Compute Fd between certain populations we're interested in
python $pyabba -f phased \
	-w $windowsize -m $minSNPs -T 2 \
	-g $inputdir"/"$inputfile \
	-o $outputdir"/fd.chr_fiji.chr_png.san_png.cla.w50m100.csv.gz" \
	-P1 $p1_2 -P2 $p2_2 -P3 $p3_2 -O $outgroup \
	--popsFile $poplabel
	
	
echo "Output file:	"$outputdir"/fd.chr_fiji.chr_png.san_png.cla.w50m100.csv.gz

Done!"