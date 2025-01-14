#!/usr/bin/bash -l

#SBATCH --job-name=vcf2bimbam
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=outscripts/%x.out
#SBATCH --error=outscripts/%x.err
#SBATCH --time=1:00:00

#### The script is used to convert annotated VCF file to bimbam file which is the required file format for ELAI analysis ####
#### Author: Wan-Ting Huang


###########################
# Arguments to pass:
workindir=${1-/scratch/sschmi13/5_elai}			
inputdir=${2-annotated_vcf/}	# Name of the directory with annotated vcf files
outputdir=${3-elai_input}	# Name of the directory with output files
popdir=${4-pop_files}	# Name of the directory with files listing out the samples of each population
chrlist=${5-chr.list}

############################################################


####################
# LOAD MODULES
####################

module load plink/1.9-beta6.27


####################
# BEGINNING OF SCRIPT
####################
cd $workindir

mkdir -p $outputdir

for chr in `cat $chrlist`
do
	echo "Working on "$chr
	
	for pop in `cat $popdir/pop.list`
	do
		plink --vcf $inputdir"/"$chr".no_lowdp_indv_chryso_sand_leuco.gatk_filtered.no_mono.biall.mac.maxmiss.mindp.maxdp.qual.annotated.vcf" \
			--recode bimbam --keep $popdir/$pop \
			--out $outputdir"/"$chr".no_lowdp_indv_chryso_sand_leuco.gatk_filtered.no_mono.biall.mac.maxmiss.mindp.maxdp.qual.${pop%%.pop}"
	done
	
	echo "Done for "$chr"
	-"
	
done
		
echo "
--------------------

Conversion to bimbam format done."	
		
		
		
