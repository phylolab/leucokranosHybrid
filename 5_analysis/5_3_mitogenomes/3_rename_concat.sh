#!/bin/bash

#SBATCH --job-name=percula_rename_concat
#SBATCH --partition=cpu
#SBATCH --output=outscripts/%x.%a.out
#SBATCH --error=outscripts/%x.%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=03:00:00

#change directory and go to mitobim output directory
cd /scratch/sschmi13/4_mitogenomes/mitobim_outputs_percula/

#initiate the loop through all samples
for i in `cat all_samples.txt`
do
cd ${i}
#keep only the last iteration directory (delete al files in the directory except the most recent one)
rm -rf `ls -t | awk 'NR>1'`
#change name of fasta entry to match the sample name
sed -i "1 s/.*/>${i}/" iteration*/${i}-percula_mtdna_assembly/${i}-percula_mtdna_d_results/${i}-percula_mtdna_out_${i}.unpadded.fasta
cd ..
done

#concatenate all samples in a single fasta file
for i in `cat all_samples.txt`
do
cat ${i}/iteration*/${i}-percula_mtdna_assembly/${i}-percula_mtdna_d_results/${i}-percula_mtdna_out_${i}.unpadded.fasta >> all_mtdna_percula.mitobim.fasta
done