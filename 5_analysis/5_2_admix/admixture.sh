#!/usr/bin/bash -l

#SBATCH --job-name=admix
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=outscripts/%x.%j.%A.out
#SBATCH --error=outscripts/%x.%j.%A.err
#SBATCH --time=1:00:00
#SBATCH --array=1-6

#### This script is used to run admixture ####
#### Author: Sarah Schmid

# WARNING /!\
# Before running this script, the input .bim file should be modified if not working with human data
# Modify .bim file because ADMIXTURE does not accept chromosome names that are not human (convert first column into 0)
# awk '{$1="0";print $0}' $file_name.bim > $file_name.bim.tmp
# mv $file_name.bim.tmp $file_name.bim 


###########################
# Arguments to pass
workindir=${1-.}
inputfile=${2-input_file/all_chr_no_lowdp_indv_chryso_sand_leuco.gatk_filtered.no_mono.biall.mac.maxmiss.mindp.maxdp.qual.bed}
outputdir=${3-admix_output}	
admix_exe=${4-./admixture/admixture}

############################################################


# Check settings

echo "Using following arguments:

Working directory:	"$workindir"

Input bed file:
	-- "$inputfile


### K value retrieved from array
k=$(sed -n ${SLURM_ARRAY_TASK_ID}p k.list)


####################
# LOAD MODULES
####################

module load htslib/1.17

####################
# BEGINNING OF SCRIPT
####################

cd $workindir

mkdir -p $outputdir

# Generate prefix name for output files
file_name=${inputfile%.bed}
prefix=${file_name##*/}

# Running admixture
$admix_exe --cv $inputfile $k > $outputdir"/log."$k".out"
 
echo "Admixture done for K="$k"!"
