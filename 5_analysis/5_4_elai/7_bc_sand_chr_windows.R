###################################|
###################################|
####
#### ancestry in windows to compare
#### to other windows statistics
#### based on ancestry inference using ELAI
#### using output generated in the cluster
#### 
#### author: Sarah Schmid
####
#### june 2022
####
####
###################################|
###################################|

# setting up environment ####

## load libraries ####
library(tidyverse)
library(ggrepel)
library(data.table)
library(GenomicRanges)
library(plyranges)
library(here)

### ancestry data ####
pop      <- read.table("../results/order_leuco_sanda.pop")
ancestry <- fread("../results/allchr_symp_mean_ps21.txt", sep = " ", h = F)
ancestry <- as_tibble(ancestry)
ancestry <- add_column(ancestry, ID = pop$V1, .before = "V1") 
ancestry <- dplyr::filter(ancestry, ID %in% c("LU009", "LU130", "LU131", "LU14", "LU21", "LU51", "LU60")) #keep only BC and F2
ancestry$ID
print("Ancestry data is loaded.")

### snps data ####
snps      <- read_delim("../results/allchr_symp_mg1_123.snpinfo.txt")  %>%
  janitor::clean_names(case = "none") %>%
  mutate(across(everything(), ~ str_replace_all(.x, " ", ""))) # to remove space in the header
snps_pos  <- snps$pos
chr       <- snps$chr
print("SNPs info is loaded")

### windows data ####
# we create a Grange object which we can use afterwards to add the ancestry info of each snps
windows_50kb <- read_csv("../results/windows_50kb.csv") %>%
  mutate(scaffold = as.numeric(str_remove(scaffold, "^Chr"))) %>%
  dplyr::rename(chr = scaffold) %>%
  makeGRangesFromDataFrame()


## prepare dataset ####

### separate ancestry values for both parents ####

#### parent 1 ancestry 
parent_1 <- ancestry %>%
  pivot_longer(cols = starts_with("V"), names_to = "SNPS_ID", values_to = "ancestry_1") %>%
  group_by(ID) %>%
  dplyr::filter(row_number() %% 2 == 0) %>%
  add_column(pos = rep(snps_pos, times = 7), .after = "SNPS_ID") %>% #times correspond to the number of individuals
  add_column(chr = rep(chr, times = 7))
print("Parent 1 dataset is prepared")

#### parent 2 ancestry
parent_2 <- ancestry %>%
  pivot_longer(cols = starts_with("V"), names_to = "SNPS_ID", values_to = "ancestry_2") %>%
  group_by(ID) %>%
  dplyr::filter(row_number() %% 2 == 1) %>%
  add_column(pos = rep(snps_pos, times = 7), .after = "SNPS_ID") %>%
  add_column(chr = rep(chr, times = 7))
print("Parent 2 dataset is prepared")

### join the two dataset and add genotype info (for the 1.8-0.2 threshold)
ancestry_final <- parent_1 %>%
  add_column(ancestry_2 = parent_2$ancestry_2) %>%
  mutate(genotype = case_when(ancestry_1 > 1.8 ~ "homo_sand",
                              ancestry_1 < 0.2 ~ "homo_chryso",
                              dplyr::between(ancestry_1, 0.8, 1.2) ~ "hetero",
                              TRUE ~ "NA")) %>%
  dplyr::rename(start = pos) %>%
  mutate(end = start, .after = start) %>%
  relocate(chr, .before = ID) %>%
  relocate(c(SNPS_ID, ID), .after = end) 

### convert to GR object ####
ancestry_GR <- makeGRangesFromDataFrame(ancestry_final, keep.extra.columns = TRUE)

### overlap SNPs and windows data for sandaracinos backcross windows ####
snp_wind <- join_overlap_inner(windows_50kb, ancestry_GR)
snp_wind_filter_san <- snp_wind %>% 
  dplyr::filter(genotype == "homo_sand") %>%
  as_tibble() %>%
  group_by(SNPS_ID) %>%
  mutate(n = dplyr::n()) %>% #number of individual with "homo_sand" per SNPs
  dplyr::filter(n == 4) %>% 
  ungroup() %>%
  select(c(seqnames, start, end)) %>%
  distinct() %>%
  write_delim("../results/allchr_BC_sand_windows.txt", delim = "\t")

### overlap SNPs and windows data for chrysopterus backcross windows ####
snp_wind_filter_chr <- snp_wind %>% 
  dplyr::filter(genotype == "homo_chryso") %>%
  as_tibble() %>%
  group_by(SNPS_ID) %>%
  mutate(n = dplyr::n()) %>% #number of individual with "homo_chr" per SNPs
  dplyr::filter(n == 1) %>% 
  ungroup() %>%
  select(c(seqnames, start, end)) %>%
  distinct() %>%
  write_delim("../results/allchr_bc_chr_windows.txt", delim = "\t")



