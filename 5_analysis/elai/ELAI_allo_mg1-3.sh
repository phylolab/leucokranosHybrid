#!/usr/bin/bash -l

#SBATCH --job-name=ELAI_allo_mg1-3
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --output=outscripts/elai_allo_mg1-3seed123.%a.out
#SBATCH --error=outscripts/elai_allo_mg1-3seed123.%a.err
#SBATCH --time=24:00:00

### The script is used to perform EALI analysis using allopatric populations as the two parental groups with admixture generation parameters as 1, 2, 3 ###

###########################
# Arguments to pass:
workindir=${1-5_ELAI}					# Working directory. Default: 5_ELAI
elai_exec=${2-ELAI/elai}
inputfileDir=${3-1_bimbam}			# Name of the directory with input files. Default: 1_bimbam
inputfileprefix=${4-_SNP.SG.Filtered.nospandel}			# Prefix of input files. Default: AllSamples_SNP.SG.Filtered.nospandel
p1=${5-chryso_fiji}
p2=${6-sanda_aus}
admixp=${7-leuco_sandadmix}
seed=${8-123}
chrlist=0_AclarkiiReference/Chromosomes.list

############################################################

############################################################
# Relaunch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $chrlist | awk '{print $1}')


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $workindir $elai_exec $inputfileDir $inputfileprefix $p1 $p2 $admixp $seed
fi

############################################################


####################
# BEGINNNING OF SCRIPT
####################
# Get the chromosome ID
chr=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${chrlist})

# Change directory
cd $workindir


# ELAI analysis
# admixture generation parameter as 1
$elai_exec -g $inputfileDir"/"$chr$inputfileprefix"."$p1".recode.geno.txt" -p 10 \
			-g $inputfileDir"/"$chr$inputfileprefix"."$p2".recode.geno.txt" -p 11 \
			-g $inputfileDir"/"$chr$inputfileprefix"."$admixp".recode.geno.txt" -p 1 \
			-pos $inputfileDir"/"$chr$inputfileprefix"."$admixp".recode.pos.txt" \
			-s 30 -C 2 -c 10 -mg 1 --ps2 -R $seed -o $chr"_allo_mg1_"$seed

# admixture generation parameter as 2
$elai_exec -g $inputfileDir"/"$chr$inputfileprefix"."$p1".recode.geno.txt" -p 10 \
			-g $inputfileDir"/"$chr$inputfileprefix"."$p2".recode.geno.txt" -p 11 \
			-g $inputfileDir"/"$chr$inputfileprefix"."$admixp".recode.geno.txt" -p 1 \
			-pos $inputfileDir"/"$chr$inputfileprefix"."$admixp".recode.pos.txt" \
			-s 30 -C 2 -c 10 -mg 2 --ps2 -R $seed -o $chr"_allo_mg2_"$seed
		
# admixture generation parameter as 3	
$elai_exec -g $inputfileDir"/"$chr$inputfileprefix"."$p1".recode.geno.txt" -p 10 \
			-g $inputfileDir"/"$chr$inputfileprefix"."$p2".recode.geno.txt" -p 11 \
			-g $inputfileDir"/"$chr$inputfileprefix"."$admixp".recode.geno.txt" -p 1 \
			-pos $inputfileDir"/"$chr$inputfileprefix"."$admixp".recode.pos.txt" \
			-s 30 -C 2 -c 10 -mg 3 --ps2 -R $seed -o $chr"_allo_mg3_"$seed

echo "---------
The script has finished."
