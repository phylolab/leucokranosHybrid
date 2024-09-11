# Introduction
This repository contains two workflows to obatin SNPs from raw sequencing data using ATLAS and GATK. We developed the workflows for [our paper](https://doi.org/10.1101/2024.03.10.584293).

# Prerequisites
## ATLAS
You can simply install the software using the below command provided by [the developers](https://anaconda.org/bioconda/atlas)
```shell
conda install bioconda::atlas
```
### For Curnagl users at the University of Lausanne
#### 1. Setting up Conda
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
#### 4. Deactivate the envirnoment
```shell
conda deactivate
```
*The scripts include a command line to activate the environment before executing ATLAS programs.*

## GATK
You can download the GATK package [here](https://github.com/broadinstitute/gatk/releases) and follow [the steps](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4) to use it.
### For Curnagl users at the University of Lausanne
Curnagl has GATK in the common software environment so we can simply put
```shell
module load openjdk/17.0.8.1_1 gatk/4.4.0.0
```
before the command lines for executing GATK programs in the bash script.

## Reference genome
The reference genome used in ATLAS workflow is [Amphiprion percula](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003047355.2/) V.1.107 and in GATK workflow is [Amphiprion clarkii](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_027123335.1/) V1. We removed the unplaced scaffolds from the downloaded files and used only chromosomes for the analysis.

# Corresponded scripts for each step
- Quality control and Trimming:
   - `Mapping_SNP_Calling_ATLAS/1_trimming.sh`
   - `Mapping_SNP_Calling_GATK/1_Run_QC_Trimming.sh`

*INDEX the reference genome using `Mapping_SNP_Calling_GATK/2A_Run_BWA_Refindexing.sh` before mapping*

- Mapping the reads:
   - `Mapping_SNP_Calling_ATLAS/2_mapping.py` for one sample
   - `Mapping_SNP_Calling_ATLAS/2_mapping_array.sh` and `Mapping_SNP_Calling_GATK/2B_Run_BWA.sh` for multiple samples
- Processing and Filtering the mapped reads:
   - `Mapping_SNP_Calling_ATLAS/3_processing.py` for one sample
   - `Mapping_SNP_Calling_ATLAS/3_processing_array.sh` and `Mapping_SNP_Calling_GATK/2C_MappingProcessing.sh` for multiple samples
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

## 2 Mapping the reads to reference genome

Using BWA version 0.7.17. The steps are :

   - Generate the Reference Indexing for BWA `bwa index`
   - Mapping `bwa mem`
   - Converting sam to bam `samtool view`
  - Sorting bam files `samtool sort`
  - Indexing sorted bam file `samtool index`
  - Statistics before filtering `bamtools stats` 

All the steps are included in the `2_mapping.py` script (adapted for python3 using Baptiste's script) and can be run using the `2_mapping_array.sh` script.

**Required input files**

* Trimmed fastq files
* Reference genome in fasta
* Reference genome index (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`) done using the `bwa index` command

## 3 Processing BAM files

The *processing BAM files* pipeline is using samtools, picard-tools and atlas and proceed to the following steps:

   - Generate Header `samtools view`
   - Add Reads groups `picard AddOrReplaceReadGroups`
   - Generate Index `samtool index`
  - Filter BAM file `samtool view -b -f 2 -F 256 -q 30`
  - Generate Index filtered bam file `samtool index`
  - Check Reads Overlap `atlas readOverlap`
  - Merge Read Groups Overlap `atlas mergeReads`
  - Check Reads Overlap after merging reads `atlas readOverlap`
  - Create Depth Mask `atlas createDepthMask`
  - Assess SoftClipping `atlas assessSoftClipping`
  - Generate Stats BAM Stats, InsertSizeMetrics `bamtools stats`
  - Generate Stats ATLAS `atlas BAMDiagnostics`and `atlas pileup`

To use atlas on Curnagl, we can use the Singularity container with the following command (be careful, in the `3_processing.py` script, we need to adapt the path of atlas)

```shell
module load singularity
export SINGULARITY_BINDPATH="/users,/work,/scratch"
singularity run /dcsrsoft/singularity/containers/atlas-0.99.sif < atlas options >
```

All the steps are included in the `3_processing.py` script and can be run using the `3_processing_array.sh` script. 

**Required input file**

* Sorted `.bam` file for each sample
* Associated `.bai` index file (no need to specify it in the command line)
