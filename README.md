## 1 Trimming reads

We trimmed the raw reads using [Trimmomatic](https://github.com/usadellab/Trimmomatic) (v.0.39, Bolger et al. 2014) and after quality control with FastQC (v.0.11.5, Andrews 2010). Adapters and low quality reads were discarded, using the following thresholds and using the adapters sequences in the TruSeq3-PE-2.fa file ([Illumina UDI adapters](https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/TruSeq/UDIndexes.htm)). 

**Required input files**

* Raw fastq files (reverse and forward if sequenced in paired-ends) for each sample

* Illumina adapters sequences (provided with Trimmomatic)

## 2 Mapping reads to reference genome

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