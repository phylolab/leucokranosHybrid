#!/usr/bin/bash -l

#SBATCH --job-name=variant_filter
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=outscripts/%x.%j.out
#SBATCH --error=outscripts/%x.%j.err
#SBATCH --time=8:00:00

#### The script is used to proceed to hard-filtering of the VCF file as suggested by GATK in absence of recalibration ####
#### as well as to additional filtering with VCFtools (optional)
#### Authors: Wan-Ting Huang, Sarah Schmid


###########################
# Arguments to pass
workindir=/scratch/sschmi13/leucokranos/4_variant_filtering					
ref_genome=/scratch/sschmi13/leucokranos/3_variant_calling/reference_genome/aclarkii_genome_chr.fna
inputdir=/scratch/sschmi13/leucokranos/3_variant_calling/3_6_joint_geno_out
inputvcf=all_chr_all_samples.vcf.gz		
outputdir=4_2_filtered_vcf
outputprefix=all_chr_all_samples

# Arguments related to VCFtools thresholds for filtering (change according to what you need)
MAC=0.02
NONMISSING=0.95
MIN_DEPTH=7
MAX_DEPTH=40

# Default setting for additional filtering
additional_filtering=true

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --additional-filtering) additional_filtering=true ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

############################################################


# Check settings
echo "Using following arguments: 

Working directory:	"$workindir"
Input directory (if relative path, from "$workindir") :	"$inputdir"
Input vcf file : "$inputvcf"
Output directory: "$outputdir"
Output prefix: "$outputprefix


echo "
If settings are wrong, please look at instructions and arguments to pass

"

####################
# LOAD MODULES
####################

module load openjdk/17.0.8.1_1 gatk/4.4.0.0 picard/3.0.0 samtools/1.17 bamtools/2.5.2 vcftools/0.1.16  bcftools/1.16

####################
# BEGINNING OF SCRIPT
####################
cd $workindir


######################
# Hard-filtering recommended by GATK in absence of VSQR (recalibration)
######################

mkdir -p $outputdir

# Check if the index file exists, if not, index the VCF file
if [[ ! -f $inputdir"/"$inputvcf".tbi" ]]; then
    echo "Indexing VCF file"
    
# Indexing VCF file
    gatk IndexFeatureFile \
         -I $inputdir"/"$inputvcf
else
    echo "Index file "$inputvcf".tbi already exists, skipping indexing step."
fi

# Extract SNPs to an individual file
echo "Extract only SNPs to a new file"

gatk SelectVariants \
	-R $ref_genome \
	-V $inputdir"/"$inputvcf \
	-select-type SNP \
	-O $outputdir"/"$inputvcf

echo "-------------"

# Annotate filter status of the SNPs
echo "Apply filters to the SNPs with these thresholds:
-- QD < 2.0
-- QUAL < 30.0
-- SOR > 3.0
-- FS > 60.0
-- MQ < 40.0
-- MQRankSum < -12.5
-- ReadPosRankSum < -8.0

Please check the thresholds!

-------------

Start to annotate the VCF file..."

gatk VariantFiltration \
    -V $outputdir"/"$inputvcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O $outputdir"/"$outputprefix".annotate.vcf.gz"
    
echo "Annotation finished!
-------------

Now we select variants passing all filters"

# Subset SNPs that pass all the filters
gatk SelectVariants \
	--exclude-filtered \
	-V $outputdir"/"$outputprefix".annotate.vcf.gz" \
	-O $outputdir"/"$outputprefix".gatk_filtered.vcf.gz"
	

echo "GATK hard-filtering done!"


######################
# Additional filtering of the VCF files with BCFtools and VCFtools (optional)
######################

if [ $additional_filtering = true ]; then

	echo "Proceeding to additional MAF, minDP, maxDP and missingness filtering..."
    
	echo "Applying following threshold for filtering: 

	MAC:	"$MAC"
	NONMISSING :	"$NONMISSING"
	MIN_DEPTH : "$MIN_DEPTH"
	MAX_DEPTH : "$MAX_DEPTH
	
	echo "
	Please check the thresholds
	
	"

    echo "Start to filter SNPs using BCFtools and VCFtools"
    echo "Keeping only biallic SNPs and removing monomorphic sites..."
    
    bcftools filter -e 'AC==0 || AC==AN' \
    $outputdir"/"$outputprefix".gatk_filtered.vcf.gz" | bcftools view -m2 -M2 -v snps -O z -o $outputdir"/"$outputprefix".gatk_filtered.no_mono.biall.vcf.gz"
    
	echo "Filtering now for depth, minor allele frequency and missingness..."
    vcftools --gzvcf $outputdir"/"$outputprefix".gatk_filtered.no_mono.biall.vcf.gz" \
            --remove-indels \
            --mac $MAC \
            --max-missing $NONMISSING \
            --minDP $MIN_DEPTH \
            --maxDP $MAX_DEPTH \
            --recode --stdout | gzip -c > $outputdir"/"$outputprefix".gatk_filtered.no_mono.biall.mac.maxmiss.mindp.maxdp.noindel.vcf.gz"

    echo "Additional VCF filtering with VCFtools done!"
fi
