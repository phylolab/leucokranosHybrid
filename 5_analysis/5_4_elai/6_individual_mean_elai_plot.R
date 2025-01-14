#### This script is used to generate karyotype plot for number of homozygous individual/sample per position/SNP ####
#### based on ancestry inference calculated by ELAI (Efficient Local Ancestry Inference; Guan 2014). ####
# To prepare the input files for this script, please follow the steps in "Local ancestry inference" of Material and Methods.

#! Can be run locally, but quite long. 

# Author: Sarah Schmid

## Load libraries --------------------------------------------------------------
library(data.table)
library(tidyverse)
library(here)


## Load data  ------------------------------------------------------------------

### Ancestry data
pop      <- read.table("../results/order_leuco_sanda.pop")
ancestry <- fread("../results/allchr_symp_mean_ps21.txt", sep = " ", h = F)
ancestry <- as_tibble(ancestry)
ancestry <- add_column(ancestry, ID = pop$V1, .before = "V1") 
ancestry$ID

### SNPs data ####
snps      <- read_delim("../results/allchr_symp_mg1_123.snpinfo.txt")  %>%
  janitor::clean_names(case = "none") %>%
  mutate(across(everything(), ~ str_replace_all(.x, " ", ""))) # to remove space in the header
snps_pos  <- snps$pos
chr       <- snps$chr
print("SNPs info is loaded")

## Prepare dataset  ------------------------------------------------------------

### Parent 1 ancestry 
parent_1 <- ancestry %>%
  pivot_longer(cols = starts_with("V"), names_to = "SNPS_ID", values_to = "ancestry_1") %>%
  group_by(ID) %>%
  dplyr::filter(row_number() %% 2 == 0) %>%
  add_column(pos = rep(snps_pos, times = 25), .after = "SNPS_ID") %>% #times correspond to the number of individuals
  add_column(chr = rep(chr, times = 25))
print("Parent 1 dataset is prepared")

### Parent 2 ancestry
parent_2 <- ancestry %>%
  pivot_longer(cols = starts_with("V"), names_to = "SNPS_ID", values_to = "ancestry_2") %>%
  group_by(ID) %>%
  dplyr::filter(row_number() %% 2 == 1) %>%
  add_column(pos = rep(snps_pos, times = 25), .after = "SNPS_ID") %>%
  add_column(chr = rep(chr, times = 25))
print("Parent 2 dataset is prepared")

### Join the two dataset and add genotype info (for the 1.8-0.2 threshold)
ancestry_final <- parent_1 %>%
  add_column(ancestry_2 = parent_2$ancestry_2) %>%
  mutate(genotype = case_when(ancestry_1 > 1.8 ~ "homo_sand",
                              ancestry_1 < 0.2 ~ "homo_chryso",
                              dplyr::between(ancestry_1, 0.8, 1.2) ~ "hetero",
                              TRUE ~ "NA"))
print("Two parent datasets are joined")


## Percentage of homozygous and heterozygous sites per individual ---------------

nb_geno_ind <- ancestry_final %>%
  ungroup() %>%
  group_by(ID, genotype) %>%
  summarise(n = n()) %>%
  mutate(percentage = (n / sum(n) * 100))  # calculate the percentage for each group


## Plotting percentage per individual ------------------------------------------

# load sample order for x-axis
pop_order <- read_csv("../data/pop_admixture_order.txt") %>% pull(1) 

# reorder the ID factord according to pop_order, so that it matches the admixture plot
nb_geno_ind <- nb_geno_ind %>%
  mutate(ID = factor(ID, levels = pop_order))


# stacked barplot
ggplot(nb_geno_ind, aes(x = ID, y = percentage, fill = genotype)) +
  geom_bar(stat = "identity") +
  scale_y_continuous() +  
  scale_fill_manual(values = c("homo_chryso" = "#942846",
                               "homo_sand" = "#47a898",  
                               "hetero" = "#dec086",     
                               "NA" = "#CCCCCC")) +     
theme_minimal() +
  labs(x = "Individual", y = "Percentage", fill = "Genotype") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank()) +
  geom_text(data = nb_geno_ind %>% 
            dplyr::filter(genotype == "homo_chryso"), 
            aes(label = scales::percent(percentage/100, accuracy = 0.1), 
            y = -1.2),  # adjust to place text below the bar
            color = "#942846",  
            size = 3, hjust = 0.5) +
  geom_text(data = nb_geno_ind %>% 
            dplyr::filter(genotype == "homo_sand"), 
            aes(label = scales::percent(percentage/100, accuracy = 0.1), 
            y = -3.2),  # adjust to place text below the bar
            color = "#47a898",  
            size = 3, hjust = 0.5)  

ggsave(file = "../plots/fig2_individual_perc.eps", width = 40, height = 30, units = "cm")
  

