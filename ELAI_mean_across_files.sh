#!/usr/bin/bash -l

#SBATCH --job-name=ELAI_mean_ps21
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=outscripts/elai_mean_ps21.out
#SBATCH --error=outscripts/elai_mean_ps21.err
#SBATCH --time=2:00:00

### The script is used to calculate the mean of estimated ancestral allele dosages for each individual at each SNP ###

###########################
# Arguments to pass:
workindir=${1-.}					# Working directory. Default: Current directory
inputfileDir=${2-5_ELAI/2_output/ps21_files}			# Name of the directory with the ps21 files. Default: 5_ELAI/2_output/ps21_files
outputfileDir=${3-5_ELAI/3_meanps21}			# Default: 5_ELAI/3_meanps21
chrlist=${4-0_AclarkiiReference/Chromosomes.list}

############################################################

####################
# BEGINNING OF SCRIPT
####################
cd $workindir

mkdir -p 5_ELAI/3_meanps21

while read -r chr
do
	if [ ! -s $outputfileDir"/"$chr"_mean_ps21.txt" ]
	then
		awk 'FNR == 1 { nfiles++; ncols = NF }
			 { for (i = 1; i <= NF; i++) sum[FNR,i] += $i
			   if (FNR > maxnr) maxnr = FNR
			 }
			 END {
				 for (line = 1; line <= maxnr; line++)
				 {
					 for (col = 1; col <= ncols; col++)
						  printf " %f", sum[line,col]/nfiles; 
					 printf "\n"
				 }
			 }' "$inputfileDir/$chr"_*.ps21.txt > $outputfileDir"/"$chr"_mean_ps21.txt"
			 
		echo "Outputted "$chr"_mean_ps21.txt
		"
	fi
done < $chrlist