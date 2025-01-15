# Introduction
This repository contains information on the pipeline, analysis and scripts used for the paper **Genomic architecture of the clownfish hybrid** ***Amphiprion leucokranos*** from Sarah Schmid, Diego A. Hartasánchez, Wan-Ting Huang, Ashton Gainsford, Geoffrey P. Jones, Nicolas Salamin. 

Authors: Sarah Schmid, Wan-Ting Huang

## Prerequisites
### Softwares
The version indicated are the one used for the manuscript analysis. However, other versions might work too. 
* Trimmomatic 0.39
* FastQC 0.12.1
* BWA 0.7.17
* Samtools 0.1.19 and 1.17   
* GATK 4.4.0.0
* Picard tools 2.9.0
* BAMtools 2.5.2
* Java 1.8.0_242
* Python 3.11.6
* bamUtil 1.0.15
* R 4.3.2
### Reference genome
The reference genome used is *[Amphiprion clarkii](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_027123335.1/)* V1. We removed the scaffolds from the downloaded files and used only the 24 full chromosomes for the analysis.
# Workflow to generate VCF files for analysis (steps 1-4)
## 1 Quality control and trimming

We trimmed the raw reads using [Trimmomatic](https://github.com/usadellab/Trimmomatic) (v.0.39, Bolger et al. 2014) and examined the quality before and after trimming using quality control tool FastQC (v.0.12.1, Andrews 2010). Adapters and low-quality reads were discarded, using the following thresholds and using the adapters sequences in the TruSeq3-PE-2.fa file ([Illumina UDI adapters](https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/TruSeq/UDIndexes.htm)). 

**Required input files**

- raw fastq files (reverse and forward if sequenced in paired-ends) for each sample
- illumina adapters sequences (provided with Trimmomatic)

**Example**

*Example for running the script on a slurm machine*

```shell
sbatch 1_run_qc_trimming.sh <samples.list> <workindir> <path_to_raw_reads> <adapters_seq.list> <min_quality> <min_length> <headcrop>
```

```shell
sbatch 1_run_qc_trimming.sh all_samples.txt . /path/to/rawdata TruSeq3-PE-2.fa 3 36 0
```

## 2 Mapping the reads to reference genome

Using [BWA-MEM](https://bio-bwa.sourceforge.net/bwa.shtml) (v0.7.17, Li et al. 2009), [SAMtools](http://www.htslib.org) (v1.17, Li et al. 2009) and [BamTools](https://hcc.unl.edu/docs/applications/app_specific/bioinformatics_tools/data_manipulation_tools/bamtools/) (v2.5.2, Barnett et al. 2011) we: 

   - indexed reference genome with `bwa index` (see `2_1_run_bwa_ref_indexing.sh` script),
   - mapped the reads with `bwa mem`,
   - converted sam to bam with `samtool view`,
   - sorted bam files with `samtool sort`,
   - indexed sorted bam file with `samtool index`,
   - generated statistics before filtering with `bamtools stats`.

Those steps are included in the the script `2_2_run_bwa_mapping.sh`. 

**Required input files**

- trimmed fastq files
- reference genome in fasta
- reference genome index (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`) done using the `bwa index` command

**Example**

Example for running the script on a slurm machine

```
sbatch 2_2_run_bwa_mapping.sh <samples.list> <workindir> <input_reads_dir> <ref_genome> <remove_sam_temp_files>
```

```shell
sbatch 2_2_run_bwa_mapping.sh all_samples.txt . 1_trimmed_reads 0_reference_genome/aclarkii_genome_chr remove_sam
```

## 3 Variant calling

### GATK preprocessing 

Before proceeding to the variant calling, we have to pre-process the mapped reads using [SAMtools](http://www.htslib.org) (v0.1.19) and [Picard Tools](http://broadinstitute.github.io/picard/) (v2.9.0) and [BamUtil](https://genome.sph.umich.edu/wiki/BamUtil) (v1.0.15). The steps are the following: 

   - generate headers with `samtools view`,
   - soft-clip beyond-end-of-reference alignments and set MAPQ to 0 for unmapped reads with `picard CleanSam`
   - add reads groups with `picard AddOrReplaceReadGroups`,
   - fix paired-reads with `samtools fixmate`
   - remove secondary alignment with `samtools view -bh -F 256` 
   - sorte the sam file with `picard SortSam` 
   - marke duplicates with `picard MarkDuplicates`
   - reorder bam files according to reference with `picard ReorderSam`
   - generate index file with `samtools index`
   - validate bam file with `bamUtil validate` 

Those steps are included in the `3_1_gatk_preprocessing.py` script which can be run on a slurm machine with the companion script `3_1_run_gatk_preprocessing.sh`. 

**Required input files**

* reference genome in fasta
* mapping files for each samples `.bam` 

**Example**

```shell
3_1_run_gatk_preprocessing.sh <samples.list> <workindir> <input_bam_dir> <ref_genome> <max_cpu> <max_mem>  
```

```shell
3_1_run_gatk_preprocessing.sh all_samples.txt . 2_bwa_mapping_out reference_genome/aclarkii_genome_chr 3 102400000
```

### Calling haplotypes

The next step is to call variants for each samples on each chromosome separatly using [GATK](https://gatk.broadinstitute.org/hc/en-us) (v4.4.0.0),  [BamTools](https://hcc.unl.edu/docs/applications/app_specific/bioinformatics_tools/data_manipulation_tools/bamtools/) (v2.5.2) and [SAMtools](http://www.htslib.org) (v0.1.19) . This is done using the script `3_3_gatk_haplotype_call.sh`. 

**Required input files**

* reference genome in fasta
* preprocessed bam files `.dedupped.reorder.bam`
* list with chromosome to process (should match the format in the reference genome)

**Example**

```shell
sbatch 3_3_gatk_haplotype_call.sh <samples.list> <workindir> <ref_genome> <input_dir_with_processed_bam> <suffix_of_processed_bam> <output_dir> <suffix_for_output_files> <list_of_chromosome_names>
```

```shell
sbatch 3_3_gatk_haplotype_call.sh all_samples.txt . reference_genome/aclarkii_genome_chr.fna 3_1_gatk_preprocessing .dedupped.reorder 3_3_gatk_haplotype .g.vcf.gz reference_genome/chromosomes.txt
```

This step generates gVCF files for each samples and each chromosome. 

### Merge gVCF files

Now, we can merge all the gVCF file generated for each chromosome together, into a single gVCF per sample using  [Picard Tools](http://broadinstitute.github.io/picard/) (v2.9.0). This is done using the script `3_4_merge_vcfs.sh`.

**Required input files**

* output from the script `3_3_gatk_haplotype_call.sh` (`.g.vcf.gz`)

**Example**

```shell
sbatch 3_4_merge_vcfs.sh <samples.list> <workindir> <input_dir_with_gvcf> <output_dir> <output_suffix>
```

```shell
sbatch 3_4_merge_vcfs.sh all_samples.txt . 3_3_gatk_haplotype 3_4_calling_output .g.vcf.gz
```

### Import gVCF into a database

The next step is to import all the gVCF files in a database, in preparation for the joint genotyping using [GATK](https://gatk.broadinstitute.org/hc/en-us) (v4.4.0.0). This is done using the script `3_5_genomic_db_import.sh`.

**Required input files**

* `g.vcf.gz` files from previous step
* reference genome in fasta
* list with chromosome to process (should match the format in the reference genome)

**Example**

```
sbatch 3_5_genomic_db_import.sh <workindir> <reference_genome> <input_gvcf_dir> <output_dir> <list_of_chromosome_names>
```

```shell
sbatch 3_5_genomic_db_import.sh . reference_genome/aclarkii_genome_chr.fna 3_4_calling_output GenomicsDB chromosomes.list
```

### Joint genotyping

We can now proceed to the joint genotyping to generate a single VCF file which includes all samples using [GATK](https://gatk.broadinstitute.org/hc/en-us) (v4.4.0.0). This is done with the script `3_6_joing_genotyping.sh`. 

**Required input files**

* reference genome in fasta
* genomic database generated with the `3_5_genomic_db_import.sh`
* list with chromosome to process (should match the format in the reference genome)

**Example**

```shell
sbatch 3_6_joint_genotyping.sh <workindir> <reference_genome> <genomic_database> <output_directory> <list_of_chromosome_names>
```

```
sbatch 3_6_joint_genotyping.sh . aclarkii_genome_chr.fna GenomicsDB 3_6_joint_geno_out chromosomes.list
```

After this step, we have a single VCF file for each chromosome. 

### Gather VCFs

We can now create a single VCF file which includes all the chromosomes using [Picard Tools](http://broadinstitute.github.io/picard/) (v2.9.0). This is done with the script `3_7_gather_vcfs.sh`. 

**Required input file**

* vcf files per chromosome generated with the `3_6_joint_genotyping.sh` script

**Example**

```shell
sbatch 3_6_joint_genotyping.sh <inputdir_with_vcfs> <output_file_name>
```

```shell
sbatch 3_6_joint_genotyping.sh 3_6_joint_geno_out all_chr_all_samples.vcf.gz
```

## 4 Variant filtering

In the absence of well-curated training resources, which are typically not available for organisms other than humans, we are unable to use VQSR (recalibration). The solution is to use hard-filtering instead. For the SNPs dataset, the following filters are recommanded by [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360037499012-I-am-unable-to-use-VQSR-recalibration-to-filter-variants)

```tex
QD < 2.0
MQ < 40.0
FS > 60.0
SOR > 3.0
MQRankSum < -12.5
ReadPosRankSum < -8.0
```

The script `4_1_variant_stats.sh` is used to generate VCF file statistics to determine what are the best threshold for filtering. 

The script `4_2_variant_filtering.sh` is used to filter the variant based on the GATK guidelines for hard-filtering, and additional filtering with VCFtools. 

# 5 Analysis

### 5.1 PCA

For the PCA, we are using plink v1.9. One of the major assumptions of PCA is that the data we use is indpendent - i.e. there are no spurious correlations among the measured variables. This is obviously not the case for most genomic data as allele frequencies are correlated due to physical linkage and linkage disequilibrium. So as a first step, we need to prune our dataset of variants that are in linkage.

Using the `--indep-pairwise` plink command, we specified a window size of 50kb and a window-step size of 10 bp (we move 10 bp each time we calculate linkage). The r2 threshold for the linkage we tolerate is 0.2, meaning that we pruned any SNPs that had a r2 above 0.2.

Check the script for the full command line `5_analysis/5_1_pca/pca_plink.sh` and the script `5_analysis/5_1_pca/pca_plotting.R` for plotting the results.

### 5.2 ADMIXTURE

We can infer indvidual ancestries with ADMIXTURE v1.3.0 (Alexander et al. 2009). To use admixture, we need a set of LD-pruned SNPs. We can use the same as for the PCA analysis. The required input is a plink `.bed` file (binary biallelic genotype table), the `.fam`and `.bim` files. The `.bim` file has to be modified because ADMIXTURE does not accept chromosome names that are not human. We will thus convert the first column into “0” using this command line

```shell
awk '{$1="0";print $0}' $file_name.bim > $file_name.bim.tmp
mv $file_name.bim.tmp $file_name.bim 
```

We can then run admixture using the script `5_2_admixture/admixture.sh`, for K=2 to K=7. ADMIXTURE produced 2 files: `.Q` which contains cluster assignments for each individual and `.P` which contains for each SNP the population allele frequencies. Then, we have to identify the best value of k clusters which is the value with lowest cross-validation error, we need to collect the cv errors.

```bash
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > $file_name.cv.error
```

The script   `5_analysis/5_2_admixture/admixture_plotting.R` is used to plot the admixture results. 

### 5.3 Mitogenomes

Scripts are in the `5_analysis/5_3_mitogenomes/` directory.

#### Merging reads

For the mitogenome reconstruction, we first have to merge overlapping paired end reads using [FLASh](https://ccb.jhu.edu/software/FLASH/)V1.2.11. This can be done using the script `1_merge_paired_reads.sh` which use has input the trimmed paired fastq files.

#### Running MITObim

Once the paired reads are merged in a single file, we can start the baiting and iterative mapping as implemented in MITOBim. We will reconstruct the mt genomes using the reference mitogenome of *A. percula* ([GenBank](https://www.ncbi.nlm.nih.gov/nuccore/KJ174497): KJ174497.1; Tao et al. 2016). To reconstruct the mitogenomes, we used the scripts `2_run_mitobim.sh`.

#### Concatenation of all fasta files

Now we can rename and concatenate all fasta files in order to then proceed to multiple alignments and to reconstruct the phylogeny. This is done using the scripts `3_rename_concat.sh`. The final concatenated fasta file is `all_mtdna_percula.mitobim.fasta`.

#### **Aligning MITObim output with MARS and MAFFT V7.505**

We can use the program [MARS](https://github.com/lorrainea/mars) (Ayad & Pissis 2017) a program which can be used in conjunction with any multiple sequence alignment program, to address the problem of circular genomes effectively and efficiently. We used the script `4_mito_alignment.sh`. 

#### **Reconstructing mtDNA tree with IQtree**

We can use [IQ-tree](http://www.iqtree.org/doc/Tutorial) v2.2.2 to reconstruct the phylogeny based on the mitochondrial genome. IQ-tree will estimate automatically select the best model and assess branch support using single branch test (SH-like approximate likelihood ratio test ([Guindon et al., 2010](https://doi.org/10.1093/sysbio/syq010)), `-alrt`) as well ultrafast bootstrap approximation (UFBoot, `-B`) ([Minh et al., 2013](https://doi.org/10.1093/molbev/mst024); [Hoang et al., 2018](https://doi.org/10.1093/molbev/msx281)). This is implemented in the script `5_run_iqtree.sh`.

#### Visualisation of the phylogeny

We can use the `.contree` file to generate a visual phylogeny with [iTOL](https://itol.embl.de/) (Letunic & Bork, 2024).

### 5.4 Local ancestry estimation with ELAI

Scripts are in the `5_analysis/5_4_elai` directory.

In order to infer genomic segments of both parental parents introgressed in *leucokranos* we used the Efficient Local Ancestry Inference (ELAI) method (Guan 2014). This method employs a two-layer Hidden Markov Model to infer the local ancestry of admixed individuals without pre-defining window sizes. It examines two layers of linkage disequilibrium—within and among the defined groups. At each variable position in the genome, it returns the most likely proportions of ancestries, where the true values are expected to be 0, 1, or 2 in a two-way admixture scenario (Seixas et al. 2018).

#### Preparing input file

We first have to prepare the input file. We have to annotate the VCF files (one VCF per chromosome) using bcftools V1.16 using the script `1_annotate_vcf.sh`. We use here the VCF consisting only in the two parental species and the hybrid (no *A. clarkii*). Then, we can convert the annotated VCF into a bimbam format using the script `2_vcf2bimbam.sh`. For each population we will have a different output file.

#### Running ELAI

We ran ELAI using the `-mg` parameters from 1 to 3. We used as parental population 1 *A. chrysopterus* from Kimbe Bay and as parental population 2 *A. sandaracinos* from Kimbe Bay. 

#### **Mean across runs and mg values**

As it was done in Seixas et al., we can calculate the means across run and different `-mg` values to have a single file for each chromosome in the end. We can use the script `4_elai_mean.sh` to do that for each chromosome on the `.ps21.txt` outputs. The final files are `ChrXX_mean_ps21.txt`. We can also generate a single file for all chromosomes using the following command

```bash
# for the allele dosage values
paste $(ls Chr*_mean_ps21.txt | sort -V) > allchr_symp_mean_ps21.txt
# for the SNPs info data
(head -n 1 $(ls Chr*_allo_mg1_123.snpinfo.txt | sort -V | head -n 1); \
 tail -n +2 -q $(ls Chr*_allo_mg1_123.snpinfo.txt | sort -V)) > allchr_symp_mg1_123.snpinfo.txt
```

#### Plotting results

We can have two different figures: one which is a summary across all individuals and one which shows the patterns for each individual.

##### Summary across all individuals

The idea here is that for each position, we count the number of individual with *sandaracinos* ancestry, the number of individual with *chyrsopterus* ancestry and the number of heterozygous. We used ggbio V1.45 to plot the karyogram. The red line on the top is whether or not there is an individual with chrysopterus ancestry at this position (we didn’t take into account the number of individuals since there is only one usually, less important then with *sandaracinos*). Here, we kept only the BC and F2 individuals (“LU009”, “LU130”, “LU131”, “LU14”, “LU21”, “LU51”, “LU60”). This is done using the script `5_plot_elai.R` which can be run on the cluster with the companion script `5_run_r_plot_elai.sh`.

##### Individual patterns

For the supplementary material, we can display the patterns for each individual *A. leucokranos*. Then, we can also calculate the percentage of homozygous *chrysopterus* windows, the percentage of homozygous *sandaracinos* windows and the percentage of heterozygous windows for each individual. This was implemented in the `6_individual_mean_elai_plot.sh` script.

##### Global patterns

For the global pattern, we counted at each position the number of individual with *sandaracinos* ancestry, the number of individual with *chyrsopterus* ancestry and the number of heterozygous. We only did that for the backcross individuals, since the full F1 individual are not informative. The visualisation was done using the R script `5_plot_elai.R`. The script can be run locally, but takes some time.

#### **Extract the backcross genomic windows**

We can extract the genomic windows in which at least 4 of the BC-sandaracinos have sandaracinos ancestry. Then, we can compare those windows with the background windows for different statistics. The script is `7_bc_sand_chr_windows.R`. The script also includes some lines to do the same for the *chrysopterus* backcross windows to not bias the “background” (if we want to compare the backcross windows with the non-backcross windows).

### 5.5 Population genomics statistics and ABBA-BABA test

All scripts are in the`5_analysis/5_5_pop_gen_stats/` directory

#### Calculating Fst, dxy and π

We used the scripts from Simon Martin ([genomics general](https://github.com/simonhmartin/genomics_general)). The script `popgenWindows.py` computes some standard population genomic statistics in sliding windows: Fst, Dxy and π. We need the `-geno` file and a `pop.file` as input. The popfile specifies the population to which each sample belongs to. We need also to specify the windows size we want to perform the calculation on. Here, we use 50kb window with a minimum number of 100 SNPs per windows. The analysis is implimented in the script `pop_gen_stats_fd.sh` and values were plotted using the scripts `plot_figure_5.R` and `plot_supp_figure_S4.R`. 

#### Calculating “barrier loci” and testing for difference between background windows and backcross/highly divergent windows

We can compare the different parental population pairs for outlier Fst values and keep only windows which are in common. This is implemented in the script `barrier_loci.R`. 

#### Are highly divergent windows biased towards some chromosomes and backcross windows identified with ELAI

##### Highly divergent windows

To detect if the number of highly divergent windows was biased toward some chromosomes, we generated a random distribution of highly divergent windows across the genome using 10,000 permutations. We then compared the observed mean value of divergent windows per chromosome to the expected mean value. We corrected the p-value for multiple testing using the Benjamini-Hochberg false-discovery rate method. We performed a two-sample Wilcoxon test to evaluate differences in FST, dxy and fdM values between the background genomic windows and the highly divergent windows. Since the statistical significance indicated by the p-value is highly influenced by the sample size and can be misleading when many observations are compared (like genomic data), we computed the effect size using Cohen’s d statistic to assess the practical significance of the difference between the two groups. We also performed a chi-square test to compare the number of barrier loci between background and divergent genomic windows. Finally, we performed a one-sample Wilcoxon rank-sum test to test whether the mean divergence value across chromosomes was statistically different from 0. This was implement in the script `test_biased_stats.R`. 

##### Backcross windows detected with ELAI

Similar to the analyses we made for the divergent windows between F1-hybrids and the parental species (see section above), we used the population genomic statistics to characterise the introgressed *A. sandaracinos* genomic tracts identified with ELAI in the *A. leucokranos* backcrosses. Here, we considered the windows with at least four out of the five *A. leucokranos* backcrosses displaying A. sandaracinos ancestry to be potential introgressed tracts. To test for a bias in some chromosomes compared to the rest of the genome, we generated a random distribution of introgressed windows across the genome using 10,000 permutations. We then compared the observed mean value of introgressed windows per chromosome to the expected mean value. We corrected the p-value for multiple testing using the Benjamini-Hochberg false-discovery rate method (Benjamini and Hochberg 1995). Then, we tested for differences in the FST, dxy and fdM values between the background and introgressed genomic windows using a two-sample Wilcoxon test and computed the effect size using Cohen’s d statistic. Finally, we performed a chi-square test to compare the number of barrier loci between background and introgressed genomic windows. This was implement in the script `test_biased_stats.R`. 

#### Mean across windows

We can calculate the mean values of each statisticts across the windows. In script `mean_across_windows.R`. 

#### Calculating fdM statistics

To check whether there are really small regions not detectable with admixture between *chrysopterus* and *sandaracinos* from PNG, we used [Simon Martin’s script](https://github.com/simonhmartin/genomics_general) `ABBABABAwindows.py`. We ran two tests, (1) one with P1 = A. san. (Christmas), P2 = A. san. (Kimbe), P3 = A. chr. (Kimbe) and the second one (2) with P1 = A. chr. (Fiji), P2 = A. chr. (Kimbe), P3 = A. san. (Kimbe). Introgression between P2 and P3 gives positive fdM values, whether introgression between P3 and P1 gives negative values. If the fdM value is different between P3 - P1 and P2 - P1, it means that there is introgression and not only ILS. By doing the two tests, we can see in which direction is the introgression (is it sandaracinos introgressing into chrysopterus, or the opposite, or in both direction similarly? We used *A. clarkii* as outgroup and used same windows size of 50kb and removed windows with less then 100 SNPs. The analysis is implimented in the script `pop_gen_stats_fd.sh`. 

##### Comparison between fdM *chrysopterus* and fdM *sandaracinos*

We can compare both fdM values for the background, q99 and q01 windows. If there is a difference, it means that introgression between the two parental species is not equal. This was done in the script `test_biased_stats.R` in the section *fdM_san vs fdM_chr in highly divergent windows*. We used a wilcoxon test and corrected for mulitple testing with bonferroni correction. We checked for effect-size using cohen’s d value. The visualisation of the results are in the script `plot_figure_4.R`.
