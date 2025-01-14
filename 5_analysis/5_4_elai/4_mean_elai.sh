#!/usr/bin/bash -l

#SBATCH --job-name=elai_mean_ps12
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=outscripts/%x.out
#SBATCH --error=outscripts/%x.err
#SBATCH --time=2:00:00

### The script is used to calculate the mean of estimated ancestral allele dosages for each individual at each SNP ###
#### Author: Wan-Ting Huang


###########################
# Arguments to pass:
workindir=${1-/scratch/sschmi13/5_elai}			
inputdir=${2-elai_output}
outputdir=${3-elai_summary_output}
chrlist=${4-chr.list}

############################################################

####################
# BEGINNING OF SCRIPT
####################
cd $workindir

mkdir -p $outputdir

while read -r chr
do
	if [ ! -s $workindir"/"$outputdir"/"$chr"_mean_ps21.txt" ]
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
			 }' "$inputdir/$chr"_*.ps21.txt > $workindir"/"$outputdir"/"$chr"_mean_ps21.txt"
			 
		echo "Outputted "$chr"_mean_ps21.txt
		"
	fi
done < $chrlist
