# scripter to generate R scripts for each chromosome in order to plot individual
# allele dosage calculated with ELAI for each A. leucokranos invididual. 
# each R script can then be run using the script 6_run_r_individual_plot.sh 

# run this script using the following command: bash 6_scripter_individual_elai_plots.sh

# author: Sarah Schmid

for i in `cat chr.list`

do echo "###################################|
###################################|
####
#### Ancestry inference using ELAI
#### Based on output generated in the cluster
#### Script to run on the cluster
####
#### June 2022
####
####
###################################|
###################################|

#### setting up environment ####
library(rlang, lib.loc = \"/users/sschmi13/R/x86_64-pc-linux-gnu-library/4.3\")
library(vctrs, lib.loc = \"/users/sschmi13/R/x86_64-pc-linux-gnu-library/4.3\")
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(data.table)
theme_set(theme_light(base_size = 12))
setwd(\"/scratch/sschmi13/5_elai/plotting_output\")

#### load data ####

#ancestry data
pop          <- read.table(\"order_leuco_sanda.pop\")
ancestry     <- fread(\"../elai_summary_output/${i}_symp_mean_ps21.txt\", h = F)

ancestry <- as_tibble(ancestry)
ancestry <- add_column(ancestry, ID = pop\$V1, .before = \"V1\")

#snps data
snps      <- read_delim(\"../elai_output/${i}_symp_mg1_123.snpinfo.txt\") %>%
	janitor::clean_names(case = \"none\") %>%
	mutate(across(everything(), ~ str_replace_all(.x, \" \", \"\"))) # to remove space in the header
snps_pos  <- snps\$pos

#parent 1 ancestry
parent_1 <- ancestry %>%
  pivot_longer(cols = starts_with(\"V\"), names_to = \"SNPS_ID\", values_to = \"ancestry_1\") %>%
  group_by(ID) %>%
  filter(row_number() %% 2 == 0) %>%
  add_column(pos = rep(snps_pos, times = 25), .after = \"SNPS_ID\")



#parent 2 ancestry
parent_2 <- ancestry %>%
  pivot_longer(cols = starts_with(\"V\"), names_to = \"SNPS_ID\", values_to = \"ancestry_2\") %>%
  group_by(ID) %>%
  filter(row_number() %% 2 == 1) %>%
  add_column(pos = rep(snps_pos, times = 25), .after = \"SNPS_ID\")

#join info on both ancestry
#for the case_when formula : for the final formula,
#we force this to evaluate as TRUE by literally using the value TRUE.
#This forces case_when to output the “ else-output-value” for any remaining
#values that weren’t previously categorized.


#with a 1.7 threshold
#ancestry_final <- list()
#for (i in 1:24) {
# ancestry_final[[i]] <- parent_1[[i]] %>%
# add_column(ancestry_2 = parent_2[[i]]\$ancestry_2) %>%
# mutate(genotype = case_when(ancestry_1 > 1.7 ~ \"homo_sand\",
#   ancestry_1 < 0.3 ~ \"homo_chryso\",
#       between(ancestry_1, 0.7, 1.3) ~ \"hetero\",
#   TRUE ~ \"NA\"))
#}

#with a 1.6 threshold
#ancestry_final <- list()
#for (i in 1:24) {
# ancestry_final[[i]] <- parent_1[[i]] %>%
# add_column(ancestry_2 = parent_2[[i]]\$ancestry_2) %>%
# mutate(genotype = case_when(ancestry_1 > 1.6 ~ \"homo_sand\",
#   ancestry_1 < 0.4 ~ \"homo_chryso\",
#       between(ancestry_1, 0.6, 1.4) ~ \"hetero\",
#   TRUE ~ \"NA\"))
#}


#with a 1.8 threshold
ancestry_final <- parent_1 %>%
  add_column(ancestry_2 = parent_2\$ancestry_2) %>%
  mutate(genotype = case_when(ancestry_1 > 1.8 ~ \"homo_sand\",
                              ancestry_1 < 0.2 ~ \"homo_chryso\",
                              between(ancestry_1, 0.8, 1.2) ~ \"hetero\",
                              TRUE ~ \"NA\"))

#with a 1.5 threshold
#ancestry_final <- list()
#for (i in 1:24) {
# ancestry_final[[i]] <- parent_1[[i]] %>%
# add_column(ancestry_2 = parent_2[[i]]\$ancestry_2) %>%
# mutate(genotype = case_when(ancestry_1 >= 1.500000 ~ \"homo_sand\",
#   ancestry_1 <= 0.500000 ~ \"homo_chryso\",
#   between(ancestry_1, 0.500001, 1.499999) ~ \"hetero\"))
#}

#colors

cols_vec <- c(\"homo_chryso\" = \"#942846\", \"homo_sand\" = \"#47A898\", \"hetero\" = \"#DEC086\", \"NA\" = \"#CCCCCC\")

#chryso: #942846
#leuco: #DEC086
#sanda: #47A898
#NA: #CCCCCC


#plotting tiles

ancestry_final %>%
    #filter(leuco_id %in% leuco_BC) %>%
    #filter(ID == \"LU009\") %>%
    ggplot(aes(x = reorder(pos, sort(as.numeric(pos))), y = ID, fill = genotype, color = genotype)) + #no white space for missing position
    #ggplot(aes(x = as.numeric(pos), y = as.factor(ID), fill = genotype, color = genotype)) +
    geom_tile(linewidth = 0.1) +
    scale_fill_manual(values = cols_vec) +
    scale_color_manual(values = cols_vec) +
    #facet_wrap(~chr, scales = \"free\") +
    theme(axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()) + 
    labs(fill = \"\", x = \"\", y = \"\") +

    ggtitle(\"${i}\")
  ggsave(paste(\"${i}\", \"_mean_ancestry_1.8.png\", sep = \"\"), height = 4, width = 12) #change here according to threshold" > "6_${i}_elai_individual_plot.R"

done