###############################|
## finding potential barrier loci in the parental species
## by looking at outlier fst windows common to the three populations
## 
## author: Sarah Schmid
##
## june 2022
################################|

# Prepare environment ####

## Load libraries ####

library(tidyverse)
library(here)

## Load data ####

fst_dxy_pi <- read_csv("../results/all_species_popgen.w50m100.csv.gz")

## Filter dataset ####

fst_dxy_pi %>%
  select(c("scaffold", "start", "end", "sites", "Fst_CHR_PNG_SAN_PNG", "Fst_CHR_FIJI_SAN_AUS", "Fst_CHR_SLM_SAN_SLM")) %>%
  mutate(q_99_png  = quantile(Fst_CHR_PNG_SAN_PNG, probs = 0.99, na.rm = T)) %>%
  mutate(q_99_slm  = quantile(Fst_CHR_SLM_SAN_SLM, probs = 0.99, na.rm = T)) %>%
  mutate(q_99_allo = quantile(Fst_CHR_FIJI_SAN_AUS, probs = 0.99, na.rm = T)) %>%
  mutate(outlier_png  =  ifelse(Fst_CHR_PNG_SAN_PNG  >= q_99_png,  "outlier", "no")) %>%
  mutate(outlier_slm  =  ifelse(Fst_CHR_SLM_SAN_SLM  >= q_99_slm,  "outlier", "no")) %>%
  mutate(outlier_allo =  ifelse(Fst_CHR_FIJI_SAN_AUS >= q_99_allo,  "outlier", "no")) %>%
  dplyr::filter(outlier_png == "outlier") %>%
  dplyr::filter(outlier_slm == "outlier") %>%
  dplyr::filter(outlier_allo == "outlier") %>%
  write_delim("../results/outlier_fst_parents.txt", delim = "\t")
  




