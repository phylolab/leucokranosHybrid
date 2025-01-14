#### This script is used to generate karyotype plot for number of homozygous individual/sample per position/SNP ####
#### based on ancestry inference calculated by ELAI (Efficient Local Ancestry Inference; Guan 2014). ####
#### To prepare the input files for this script, please follow the steps in "Local ancestry inference" of Material and Methods.

#### Author: Sarah Schmid


#! Can be run locally, but quite long. 


## Load libraries
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(data.table)
library(ggbio)
library(GenomicRanges)
library(grDevices)
library(here)


## Load data
### A.clarkii karyotype 
clarkii.karyotype      <- read.table("../results/karyotype_clarkii.txt", h=T)
clarkii.GR             <- makeGRangesFromDataFrame(clarkii.karyotype) # Convert to Granges
seqlengths(clarkii.GR) <- clarkii.karyotype$End # Add sequence length info
print("Setting up karyotype done!")

### Ancestry data ####
pop      <- read.table("../results/order_leuco_sanda.pop")
ancestry <- fread("../results/allchr_symp_mean_ps21.txt", sep = " ", h = F)
ancestry <- as_tibble(ancestry)
ancestry <- add_column(ancestry, ID = pop$V1, .before = "V1") 
ancestry <- dplyr::filter(ancestry, ID %in% c("LU009", "LU130", "LU131", "LU14", "LU21", "LU51", "LU60")) #keep only BC and F2
ancestry$ID
print("Ancestry data is loaded.")

### SNPs data ####
snps      <- read_delim("../results/allchr_symp_mg1_123.snpinfo.txt")  %>%
  janitor::clean_names(case = "none") %>%
  mutate(across(everything(), ~ str_replace_all(.x, " ", ""))) # to remove space in the header
snps_pos  <- snps$pos
chr       <- snps$chr
print("SNPs info is loaded")


## Prepare dataset
### Separate ancestry values for both parents

#### Parent 1 ancestry 
parent_1 <- ancestry %>%
  pivot_longer(cols = starts_with("V"), names_to = "SNPS_ID", values_to = "ancestry_1") %>%
  group_by(ID) %>%
  dplyr::filter(row_number() %% 2 == 0) %>%
  add_column(pos = rep(snps_pos, times = 7), .after = "SNPS_ID") %>% #times correspond to the number of individuals
  add_column(chr = rep(chr, times = 7))
print("Parent 1 dataset is prepared")

#### Parent 2 ancestry
parent_2 <- ancestry %>%
  pivot_longer(cols = starts_with("V"), names_to = "SNPS_ID", values_to = "ancestry_2") %>%
  group_by(ID) %>%
  dplyr::filter(row_number() %% 2 == 1) %>%
  add_column(pos = rep(snps_pos, times = 7), .after = "SNPS_ID") %>%
  add_column(chr = rep(chr, times = 7))
print("Parent 2 dataset is prepared")

### Join the two dataset and add genotype info (for the 1.8-0.2 threshold)
ancestry_final <- parent_1 %>%
  add_column(ancestry_2 = parent_2$ancestry_2) %>%
  mutate(genotype = case_when(ancestry_1 > 1.8 ~ "homo_sand",
                              ancestry_1 < 0.2 ~ "homo_chryso",
                              dplyr::between(ancestry_1, 0.8, 1.2) ~ "hetero",
                              TRUE ~ "NA"))
print("Two parent datasets are joined")

### Number of genotype type per position (for one SNP, how  many individuals homo_san, homo_chr, NA and hetero)
nb_geno <- ancestry_final %>%
  ungroup() %>%
  group_by(chr, pos, genotype, .drop = FALSE) %>%
  dplyr::summarise(n = n()) %>%
  pivot_wider(names_from = genotype, values_from = n) %>%
  ungroup() %>%
  dplyr::rename(start = pos) %>%
  mutate(end = start, .after = start) %>%
  select(chr, start, end, hetero, homo_sand, homo_chryso) %>%
  mutate(homo_chryso = replace_na(homo_chryso, 0)) %>%
  mutate(homo_sand = replace_na(homo_sand, 0))
print("Finish counting the genotype")

### Convert to GRanges object
nb_geno_GR             <- makeGRangesFromDataFrame(nb_geno, keep.extra.columns = TRUE)
seqlengths(nb_geno_GR) <- seqlengths(clarkii.GR)[names(seqlengths(nb_geno_GR))]

### Granges object with only position which are not 0 for homo_chryso
geno_chryso <- nb_geno_GR[nb_geno_GR$homo_chryso != '0']
print("Generate a Granges object with positions that are homozygous for A.chrysopterus")

## Plotting karyotype
png(file="../plots/allchr_symp_ELAI_geno_count.png", width=5000, height=3500)
ggbio::autoplot(nb_geno_GR, layout = "karyogram", aes(color = homo_sand, fill = homo_sand)) +
  scale_color_gradient(low = "#ffffff", high = "#1f594e") +
  scale_fill_gradient(low = "#ffffff", high = "#1f594e") +
  layout_karyogram(data = geno_chryso, geom = "rect", ylim = c(11, 10.3), color = "#942846") +
  theme_null()
dev.off()

print("The script has finished!")




