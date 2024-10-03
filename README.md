# Introduction
This repository contains the scripts used to obtain SNP data from raw sequencing data as well as scripts used for subsequent analysis and figures for the following publication [our paper](https://doi.org/10.1101/2024.03.10.584293).
# Prerequisites - ATLAS
You can simply install the software using the below command provided by [the developers](https://anaconda.org/bioconda/atlas) or follow [the steps](https://atlaswiki.netlify.app/getting_started) to install from scratch.
```shell
conda install bioconda::atlas
```
### For Curnagl users at the University of Lausanne
#### 1. Setting up conda
```shell
module load miniconda3
conda init bash
```
`conda init bash` will hang on a sudo password input, just ignore it (ctrl-c). For more details, please consult this [wiki page](https://wiki.unil.ch/ci/books/high-performance-computing-hpc/page/using-conda-and-anaconda).
#### 2. Create an environment
```shell
conda create --ATLAS
```
#### 3. Install the dependencies in this environment
```shell
conda activate ATLAS
conda install -c conda-forge -c bioconda atlas
```
#### 4. Deactivate the environment
```shell
conda deactivate
```
*The scripts include a command line to activate the environment before executing ATLAS programs.*
### Reference genome
The reference genome used in ATLAS workflow is [Amphiprion percula](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003047355.2/) V.1.107 . 
## Prerequisites - GATK pipeline
### Softwares
The version indicated are the one used for the manuscript analysis. However, other versions might work too. 
* Trimmomatic 0.39
* FastQC 0.12.1
* BWA 0.7.17
* Samtools 1.17   
* GATK 4.4.0.0
* Picard tools 2.9.0
* BAMtools 2.5.2
* Java 1.8.0_242
* Python 3.11.6
* bamUtil 1.0.15
* R 4.3.2
### Reference genome
The reference genome used is *[Amphiprion clarkii](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_027123335.1/)* V1. We removed the unplaced scaffolds from the downloaded files and used only chromosomes for the analysis.
# Summary of corresponding scripts for each step
**Quality control and trimming of reads with Trimmomatic V0.39
   - `atlas_pipeline/1_trimming.sh`
   - `1_filtering_reads/1_run_qc_trimming.sh`
**Indexing reference genome before mapping with BWA V0.7.17**
* `2_mapping/2_1_run_bwa_ref_indexing.sh`
**Mapping the reads with BWA V0.7.17**
   - `atlas_pipeline/2_mapping.py` for one sample
   - `atlas_pipeline/2_mapping_array.sh` for multiple samples 
   - `2_mapping/2_2_run_bwa_mapping.sh` 
**GATK preprocessing of mapped reads**
   - `atlas_pipeline/3_processing.py` for one sample
   - `atlas_pipeline/3_processing_array.sh` for multiple samples
   - `3_variant_calling/3_1_gatk_preprocessing.py` for one sample and use the script `3_1_run_gatk_preprocessing.sh` to run the script on Slurm for multiple samples. 
**SNPs calling**
* *
- Calling SNPs:
   - `Mapping_SNP_Calling_GATK/3*`
- Filtering SNPs:
   - `Mapping_SNP_Calling_GATK/4*`

# Workflow in details
## 1. Quality control and Trimming

We trimmed the raw reads using [Trimmomatic](https://github.com/usadellab/Trimmomatic) (v.0.39, Bolger et al. 2014) and examined the quality before and after trimming using quality control tool FastQC (v.0.12.1, Andrews 2010). Adapters and low-quality reads were discarded, using the following thresholds and using the adapters sequences in the TruSeq3-PE-2.fa file ([Illumina UDI adapters](https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/TruSeq/UDIndexes.htm)). 

**Required input files**

- Raw fastq files (reverse and forward if sequenced in paired-ends) for each sample
- Illumina adapters sequences (provided with Trimmomatic)

**Example**

*Here use `1_Run_QC_Trimming.sh`*

```shell
sbatch 1_Run_QC_Trimming.sh AllSamples.txt . /Path/To/RawData Mapping_SNP_Calling_GATK/TruSeq3-PE-2.fa 20 50 10
```

## 2. Mapping the reads to reference genome

Using BWA-MEM (v0.7.17, Li et al. 2009), SAMtools (v1.17, Li et al. 2009) and BamTools (v2.5.2, Barnett et al. 2011), we 

   - indexed reference genome by `bwa index`,
   - mapped the reads by `bwa mem`,
   - converted sam to bam by `samtool view`,
   - sorted bam files by `samtool sort`,
   - indexed sorted bam file by `samtool index`,
   - generated statistics before filtering by `bamtools stats`.

These steps are included in `2_mapping.py`, `2_mapping_array.sh`, and `2B_Run_BWA.sh`.

**Required input files**

- Trimmed fastq files
- Reference genome in fasta
- Reference genome index (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`) done using the `bwa index` command

**Example**

*Here use `Mapping_SNP_Calling_GATK/2B_Run_BWA.sh`*

```shell
sbatch 2A_Run_BWA_RefIndexing.sh 
sbatch 2B_Run_BWA.sh AllSamples.txt . 1_TrimmedReads 0_AclarkiiReference/AclarkiiGenome.Chr RemoveSAM
```

## 3. Processing and Filtering the mapped reads

Using SAMtools, [Picard Tools](http://broadinstitute.github.io/picard/)(v3.0.0) and ATLAS (v0.9, Link et al. 2017), we

   - generated headers by `samtools view`,
   - added reads groups by `picard AddOrReplaceReadGroups`,
   - generated index by `samtool index`,
   - filtered BAM files by `samtool view -b -f 2 -F 256 -q 30`,
   - generated index of filtered BAM files by `samtool index`,
   - checked overlapping reads by `atlas readOverlap`,
   - merged reads with removing overlapping reads by `atlas mergeReads`,
   - checked overlapping reads after merging by `atlas readOverlap`,
   - created depth mask by `atlas createDepthMask`,
   - assessed SoftClipping by `atlas assessSoftClipping`,
   - generated statistics of merged reads by `picard CollectInsertSizeMetrics`, `bamtools stats`, `atlas BAMDiagnostics` and `atlas pileup`.

These steps are included in `3_processing.py`, `3_processing_array.sh`, and `2C_MappingProcessing.sh`.


**Required input files**

- Sorted `.bam` file for each sample
- Associated `.bai` index file (no need to specify it in the command line)

**Example**

*Here use `Mapping_SNP_Calling_GATK/2C_MappingProcessing.sh`*

```shell
sbatch 2C_MappingProcessing.sh AllSamples.txt . 2B_MappingOutput 30 atlas
```


To use atlas on Curnagl, we can use the Singularity container with the following command (be careful, in the `3_processing.py` script, we need to adapt the path of atlas)

```shell
module load singularity
export SINGULARITY_BINDPATH="/users,/work,/scratch"
singularity run /dcsrsoft/singularity/containers/atlas-0.99.sif < atlas options >
```


## 4. Calling SNPs
### GATK
Following GATK Best Practices recommendations (DePristo et al., 2011; Van der Auwera & O'Connor 2020) to call the variants, we 

   - indexed reference genome by `gatk CreateSequenceDictionary` and `samtools faidx`, creating reuquired `.dict` and `.fai` files,
   - called variants chromosome by chromosome by `gatk HaplotypeCaller`,
   - merged variants of chromosomes for each sample by `picard MergeVcfs`,
   - imported all merged gVCF files into a database workspace by `gatk GenomicsDBImport`,
   - performed joint-calling by `gatk GenotypeGVCFs`.

These steps are executed chronologically by scripts 3A to 3E.

**Required input files**

- Reference genome in fasta
- Merged `.bam` file for each sample
- Associated `.bai` index file (no need to specify it in the command line)

**Example**

```bash
sbatch 3A_CreateIndexFilesForGATK.sh
sbatch 3B_HaplotypeCaller_ByChromosomes.sh AllSamples.txt . 0_AclarkiiReference/AclarkiiGenome.Chr.fna 2C_FilteredMapping .BWA.Aclarkii.Sort.Filt_mergedReads 3B_CallingHaplotype g.vcf.gz 0_AclarkiiReference/Chromosomes.txt
sbatch 3C_MergeVCFs.sh AllSamples.txt . 0_AclarkiiReference/AclarkiiGenome.Chr.fna 3B_CallingHaplotype 3C_CallingOutput .g.vcf.gz
sbatch 3D_GenomicsDBImport.sh
sbatch 3E_JointGenotyping.sh
```

## 5. Filtering SNPs
Please use `4A_VariantStatistics.sh` to calculate the statistics of your variants, which can be visualised by `VariantStats.R` and will help you determine the thresholds. Next you can use `4B_Filtering.sh` to filter the varaints and check the statistics again. Using VCFtools (v0.1.16, Danecek et al. 2011), we

   - filtered the variants by `vcftools`,
   - piped the standard output of `vcftools` with `gzip -c`,
   - saved explicitly as a new file by `>`.


**Required input files**

- Files with variants


**Example**

```bash
sbatch 4A_VariantStatistics.sh
sbatch 4B_Filtering.sh
```
