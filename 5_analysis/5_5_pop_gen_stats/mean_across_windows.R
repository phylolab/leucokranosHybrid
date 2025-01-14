###############################|
## 
## mean values of fst, dxy, π and TD 
## in 50kb genomic windows
##
## author: Sarah Schmid
##
## june 2022
################################|

# prepare environment ####

## load libraries ####

library(ggrepel)
library(tidyverse)
library(patchwork)
library(here)

## load data ####
fst_dxy_pi    <- read_csv("../results/all_species_popgen.w50m100.csv.gz")

# calculate statistics across genome ####

## mean across windows ####
mean_pi <- fst_dxy_pi %>% 
  select(starts_with("pi")) %>%
  dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = starts_with("pi"), names_to = "pop", values_to = "mean_pi", names_prefix = "pi_")

mean_dxy <- fst_dxy_pi %>% 
  select(starts_with("dxy")) %>%
  dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = starts_with("dxy"), names_to = "pop", values_to = "mean_dxy", names_prefix = "dxy_")

mean_Fst <- fst_dxy_pi %>% 
  select(starts_with("Fst")) %>%
  dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = starts_with("Fst"), names_to = "pop", values_to = "mean_Fst", names_prefix = "Fst_")


## sd across windows ####
sd_pi <- fst_dxy_pi %>% 
  select(starts_with("pi")) %>%
  dplyr::summarise(across(where(is.numeric), ~ sd(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = starts_with("pi"), names_to = "pop", values_to = "sd_pi", names_prefix = "pi_")

sd_dxy <- fst_dxy_pi %>% 
  select(starts_with("dxy")) %>%
  dplyr::summarise(across(where(is.numeric), ~ sd(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = starts_with("dxy"), names_to = "pop", values_to = "sd_dxy", names_prefix = "dxy_")

sd_Fst <- fst_dxy_pi %>% 
  select(starts_with("Fst")) %>%
  dplyr::summarise(across(where(is.numeric), ~ sd(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = starts_with("Fst"), names_to = "pop", values_to = "sd_Fst", names_prefix = "Fst_")


## new table and export it to csv ####
mean_sd_fst_dxy <- left_join(mean_Fst, mean_dxy, by = "pop") %>%
  left_join(sd_Fst, by = "pop") %>%
  left_join(sd_dxy, by = "pop") %>%
  write_delim("../results/mean_values/mean_sd_fst_dxy_50kb.txt", delim = "\t")

mean_sd_pi <- left_join(mean_pi, sd_pi, by = "pop") %>%
  write_delim("../results/mean_values/mean_sd_pi.txt", delim = "\t")


## extract info for supp. mat. ##
fst_dxy_pi %>%
  select(!starts_with(c("pi", "dxy"))) %>%
  select(!c("mid")) %>%
  select(!contains("miss")) %>%
  write_delim("../results/windows_values/fst_windows.txt", delim = "\t")

fst_dxy_pi %>%
  select(!starts_with(c("pi", "fst"))) %>%
  select(!c("mid")) %>%
  select(!contains("miss")) %>%
  write_delim("../results/windows_values/dxy_windows.txt", delim = "\t")

fst_dxy_pi %>%
  select(!starts_with(c("fst", "dxy"))) %>%
  select(!c("mid")) %>%
  select(!contains("miss")) %>%
  write_delim("../results/windows_values/pi_windows.txt", delim = "\t")


