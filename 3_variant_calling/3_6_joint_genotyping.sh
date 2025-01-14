#!/bin/bash

#SBATCH --job-name=JointGenotyping
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --output=outscripts/JointGenotyping_chr%a.out
#SBATCH --error=outscripts/JointGenotyping_chr%a.err
#SBATCH --time=24:00:00

#### The script is used to perform GATK GenotypeGVCFs to create a single VCF file per chromosome with all samples included ####
#### Author: Wan-Ting Huang




############################################################
# Arguments to pass
workindir=${1-/scratch/sschmi13/leucokranos/2_bwa_mapping}					
ref_genome=${2-reference_genome/aclarkii_genome_chr.fna}
db=${3-GenomicsDB} # name of database create with 3_5_genomics_db_import.sh. Should be in workindir
outputfile_dir=${4-3_6_joint_geno_out}	 	
file_chromosomes=${5-reference_genome/chromosomes.list}  # file with chromosomes names


############################################################
## relaunch the script as array with the number of jobs corresponding to the number of files
nline=$(wc -l $file_chromosomes | awk '{print $1}')


if [[ "$SLURM_ARRAY_TASK_ID" == "" ]]; then
     # Relaunch this script as an array
     exec sbatch --array=1-$nline $0 $workindir $ref_genome $db $outputfile_dir $file_chromosomes 
fi


############################################################
# check settings
# get name of samples for array
chrom=$(sed -n ${SLURM_ARRAY_TASK_ID}p $file_chromosomes)

echo "Using following arguments: 

Working directory:	"$workindir"
Name of database :	"$db"
Output directory: "$outputfile_dir


echo "
If settings are wrong, please look at instructions and arguments to pass

"

####################
# LOAD MODULES
####################

module load openjdk/17.0.8.1_1 gatk/4.4.0.0 picard/3.0.0 samtools/1.17 bamtools/2.5.2 

####################
# BEGINNING OF SCRIPT
####################

cd $workindir

# perform joint genotyping on GenomicsDB workspace 
echo "Start joint genotyping for "$chrom"..."

gatk --java-options "-Xmx4g -XX:ParallelGCThreads=2" GenotypeGVCFs \
	-R $workindir"/"$ref_genome \
	-V gendb://$db \
	-O $workindir"/"$outputfile_dir"/"$chrom"_all_samples.vcf.gz" \
	-L $chrom
	
echo "Joint genotyping is done and output a file "$outputfile_dir"/all_samples.vcf.gz."

####################
# END
####################

echo "Joint genotyping for "$chrom "is done! You can now run the script 3_7_gathervcfs.sh to merge all chromosomes together."


