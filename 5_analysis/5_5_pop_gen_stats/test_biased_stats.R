###############################|
## fst and dxy value shifted towards sandaracinos in F1?
## fst and dxy outlier values more present in some chromsomes?
## A. leucokranos
##
## author: Sarah Schmid
##
## june 2022
################################|

# importance of effect size: https://www.scribbr.com/statistics/effect-size
# wilcoxon-test: https://www.datanovia.com/en/lessons/wilcoxon-test-in-r/


## Load libraries ####

library(tidyverse)
library(rstatix)
library(patchwork)
library(coin)
library(janitor)
library(multcompView)
library(here)

## Environment ####

theme_set(theme_light(base_size = 12))
options(digits = 2)
color.list <- c("#B33771", "#FC427B", "#BDC581", "#25CCF7", "#F97F51", "#6D214F", "#3B3B98", "#EAB543")


## Load data ####
fst_dxy_pi    <- read_csv("../results/all_species_popgen.w50m100.csv.gz") %>% mutate(scaffold = as.numeric(str_remove(scaffold, "^Chr")))
fdm_chr       <- read_csv("../results/fd.chr_fiji.chr_png.san_png.cla.w50m100.csv.gz") %>% 
  mutate(scaffold = as.numeric(str_remove(scaffold, "^Chr"))) %>%
  dplyr::rename(fdM_chr = fdM) %>%
  select(scaffold, start, end, fdM_chr)
fdm_san       <- read_csv("../results/fd.san_aus.san_png.chr_png.cla.w50m100.csv.gz") %>%
  mutate(scaffold = as.numeric(str_remove(scaffold, "^Chr"))) %>%
  dplyr::rename(fdM_san = fdM) %>%
  select(scaffold, start, end, fdM_san)


bc_sand       <- read_delim("../../5_5_elai/results/allchr_bc_sand_windows.txt") %>%
  dplyr::rename("scaffold" = seqnames) %>%
  mutate(bc_sand_window = rep("bc_sand_win"))

bc_chryso     <- read_delim("../../5_5_elai/results/allchr_bc_chr_windows.txt") %>%
  dplyr::rename("scaffold" = seqnames) %>%
  mutate(bc_chryso_window = rep("bc_chryso_win"))


## Combine data in a single dataset ####
df_all_stats <- left_join(fst_dxy_pi, fdm_san) %>%
  left_join(fdm_chr) %>%
  left_join(bc_sand) %>%
  mutate(bc_sand_window = replace_na(bc_sand_window, "background")) %>%
  left_join(bc_chryso) %>%
  mutate(bc_chryso_window = replace_na(bc_chryso_window, "background"))


### All windows dataset ####
# We want a dataset with all the windows we analysed, that we will use as universe 
all_windows <- data.frame(fst_dxy_pi) %>% select(c(scaffold, start, end))

## Compute outlier windows for fst, dxy, fdM ####
df_all_outlier <- df_all_stats %>%
  drop_na() %>%
  #outlier for fst (barrier loci, highly differentiated SNPs among parental populations)
  mutate(q_99_fst_png  = quantile(Fst_CHR_PNG_SAN_PNG, probs = 0.99, na.rm = T)) %>%
  mutate(q_99_fst_slm  = quantile(Fst_CHR_SLM_SAN_SLM, probs = 0.99, na.rm = T)) %>%
  mutate(q_99_fst_allo = quantile(Fst_CHR_FIJI_SAN_AUS, probs = 0.99, na.rm = T)) %>%
  mutate(outlier_fst_png  =  ifelse(Fst_CHR_PNG_SAN_PNG  >= q_99_fst_png,  "outlier", "no")) %>%
  mutate(outlier_fst_slm  =  ifelse(Fst_CHR_SLM_SAN_SLM  >= q_99_fst_slm,  "outlier", "no")) %>%
  mutate(outlier_fst_allo =  ifelse(Fst_CHR_FIJI_SAN_AUS >= q_99_fst_allo,  "outlier", "no")) %>%
  mutate(barrier_loci     =  ifelse(outlier_fst_png ==  "outlier" & outlier_fst_slm == "outlier" & outlier_fst_allo == "outlier", "barrier", "no")) %>%
  #outlier for fdm san
  mutate(q_99_fdM_san = quantile(fdM_san, probs = 0.99, na.rm = T)) %>%
  mutate(q_01_fdM_san = quantile(fdM_san, probs = 0.01, na.rm = T)) %>%
  mutate(outlier_q99_fdM_san  =  ifelse(fdM_san  >= q_99_fdM_san,  "outlier", "no")) %>%
  mutate(outlier_q01_fdM_san  =  ifelse(fdM_san <= q_01_fdM_san,  "outlier", "no")) %>%
  #outlier for fdm chr
  mutate(q_99_fdM_chr = quantile(fdM_chr, probs = 0.99, na.rm = T)) %>%
  mutate(q_01_fdM_chr = quantile(fdM_chr, probs = 0.01, na.rm = T)) %>%
  mutate(outlier_q99_fdM_chr =  ifelse(fdM_chr  >= q_99_fdM_chr,  "outlier", "no")) %>%
  mutate(outlier_q01_fdM_chr =  ifelse(fdM_chr  <= q_01_fdM_chr,  "outlier", "no"))

  
## Compute difference #####

### Difference in Fst #####
df_final <- df_all_outlier %>%
  drop_na() %>%
  #calculate difference
  mutate(diff_bc_png = Fst_CHR_PNG_leuco_BC - Fst_leuco_BC_SAN_PNG) %>%
  mutate(diff_bc_slm = Fst_CHR_SLM_leuco_BC - Fst_leuco_BC_SAN_SLM) %>% 
  mutate(diff_bc_allo = Fst_CHR_FIJI_leuco_BC - Fst_leuco_BC_SAN_AUS) %>% 
  mutate(diff_f1_png  = Fst_CHR_PNG_leuco_F1 - Fst_leuco_F1_SAN_PNG) %>%
  mutate(diff_f1_slm  = Fst_CHR_SLM_leuco_F1 - Fst_leuco_F1_SAN_SLM) %>% 
  mutate(diff_f1_allo = Fst_CHR_FIJI_leuco_F1 - Fst_leuco_F1_SAN_AUS) %>% 
  pivot_longer(cols = starts_with("diff"), names_to = "pop", values_to = "fst_diff", names_prefix = "diff_") %>%
  #calculate quantile
  group_by(pop) %>%
  mutate(q99 = quantile(fst_diff, probs = 0.99, na.rm = T)) %>% 
  mutate(q01 = quantile(fst_diff, probs = 0.01, na.rm = T)) %>%
  #define outlier windows
  mutate(outlier_q99  = ifelse(fst_diff >= q99, "outlier", "no")) %>%
  mutate(outlier_q01  = ifelse(fst_diff <= q01, "outlier", "no")) %>%
  mutate(BG_or_out_diff_F1 = case_when(fst_diff >= q99 ~ "out_q99",
                                       fst_diff <= q01 ~ "out_q01",
                                       TRUE ~ "background")) %>%
  #global median and mean fst difference 
  mutate(median_glob_fst_diff = median(fst_diff, na.rm = T)) %>%
  mutate(mean_glob_fst_diff = mean(fst_diff, na.rm = T)) %>%
  #median fst difference per chromosome
  ungroup() %>%
  group_by(pop, scaffold) %>%
  mutate(median_chr_fst_diff = median(fst_diff, na.rm = T)) %>%
  mutate(mean_chr_fst_diff = mean(fst_diff, na.rm = T)) %>%
  ungroup()


### Difference in Dxy #####
df_final <- df_all_outlier %>%
  drop_na() %>%
  #calculate difference
  mutate(diff_bc_png = dxy_CHR_PNG_leuco_BC - dxy_leuco_BC_SAN_PNG) %>%
  mutate(diff_bc_slm = dxy_CHR_SLM_leuco_BC - dxy_leuco_BC_SAN_SLM) %>% 
  mutate(diff_bc_allo = dxy_CHR_FIJI_leuco_BC - dxy_leuco_BC_SAN_AUS) %>% 
  mutate(diff_f1_png  = dxy_CHR_PNG_leuco_F1 - dxy_leuco_F1_SAN_PNG) %>%
  mutate(diff_f1_slm  = dxy_CHR_SLM_leuco_F1 - dxy_leuco_F1_SAN_SLM) %>% 
  mutate(diff_f1_allo = dxy_CHR_FIJI_leuco_F1 - dxy_leuco_F1_SAN_AUS) %>% 
  pivot_longer(cols = starts_with("diff"), names_to = "pop", values_to = "dxy_diff", names_prefix = "diff_") %>%
  #calculate quantile
  group_by(pop) %>%
  mutate(q99 = quantile(dxy_diff, probs = 0.99, na.rm = T)) %>% 
  mutate(q01 = quantile(dxy_diff, probs = 0.01, na.rm = T)) %>%
  #define outlier windows
  mutate(outlier_q99  = ifelse(dxy_diff >= q99, "outlier", "no")) %>%
  mutate(outlier_q01  = ifelse(dxy_diff <= q01, "outlier", "no")) %>%
  mutate(BG_or_out_diff_F1 = case_when(dxy_diff >= q99 ~ "out_q99",
                                       dxy_diff <= q01 ~ "out_q01",
                                       TRUE ~ "background")) %>%
  #global median and mean dxy difference 
  mutate(median_glob_dxy_diff = median(dxy_diff, na.rm = T)) %>%
  mutate(mean_glob_dxy_diff = mean(dxy_diff, na.rm = T)) %>%
  #median dxy difference per chromosome
  ungroup() %>%
  group_by(pop, scaffold) %>%
  mutate(median_chr_dxy_diff = median(dxy_diff, na.rm = T)) %>%
  mutate(mean_chr_dxy_diff = mean(dxy_diff, na.rm = T)) %>%
  ungroup()

#########################################################################################|
# comparison between backcross windows/outlier diff. F1-parents and background windows ####
#########################################################################################|

## fdM, Fst, dxy and pi in outlier diff fst F1-parents vs background windows ####

### prepare df ####
df_fst <- df_final %>%
  select(c(scaffold, start, end, outlier_q99, outlier_q01, BG_or_out_diff_F1, starts_with("Fst"))) %>%
  pivot_longer(cols = starts_with("Fst"), names_to = "pop_Fst", values_to = "Fst", names_prefix = "Fst_")

df_dxy <- df_final %>%
  select(c(scaffold, start, end, outlier_q99, outlier_q01, BG_or_out_diff_F1, starts_with("dxy"))) %>%
  pivot_longer(cols = starts_with("dxy"), names_to = "pop_dxy", values_to = "dxy", names_prefix = "dxy_")

df_pi <- df_final %>%
  select(c(scaffold, start, end, outlier_q99, outlier_q01, BG_or_out_diff_F1, starts_with("pi"))) %>%
  pivot_longer(cols = starts_with("pi"), names_to = "pop_pi", values_to = "pi", names_prefix = "pi_")

df_fdM_san <- df_final %>%
  select(c(scaffold, start, end, outlier_q99, outlier_q01, BG_or_out_diff_F1, starts_with("fdM_san"))) %>%
  distinct(scaffold, start, end, outlier_q99, outlier_q01, fdM_san, BG_or_out_diff_F1)

df_fdM_chr <- df_final %>%
  select(c(scaffold, start, end, outlier_q99, outlier_q01, BG_or_out_diff_F1, starts_with("fdM_chr"))) %>%
  distinct(scaffold, start, end, outlier_q99, outlier_q01, fdM_chr, BG_or_out_diff_F1)

# for comparison between the two fdM values
df_fdM_san_chr <- df_final %>%
  select(c(scaffold, start, end, outlier_q99, outlier_q01, BG_or_out_diff_F1, starts_with("fdM"))) %>%
  distinct(scaffold, start, end, outlier_q99, outlier_q01, fdM_chr, fdM_san, BG_or_out_diff_F1)

### extract mean values ####
df_fst %>% group_by(outlier_q99, outlier_q01, pop_Fst) %>%
  dplyr::summarise(mean_Fst = mean(Fst), sd_Fst = sd(Fst)) %>%
  write_delim("../results/mean_highly_div_windows/mean_Fst_diff_fst_f1_parent_windows.txt", delim = "\t")

df_dxy %>% group_by(outlier_q99, outlier_q01, pop_dxy) %>%
  dplyr::summarise(mean_dxy = mean(dxy), sd_dxy = sd(dxy)) %>%
  write_delim("../results/mean_highly_div_windows/mean_dxy_diff_fst_f1_parent_windows.txt", delim = "\t")

df_pi %>% group_by(outlier_q99, outlier_q01, pop_pi) %>%
  dplyr::summarise(mean_pi = mean(pi), sd_pi = sd(pi)) %>%
  write_delim("../results/mean_highly_div_windows/mean_pi_diff_fst_f1_parent_windows.txt", delim = "\t")

df_fdM_san %>% group_by(outlier_q99, outlier_q01) %>%
  dplyr::summarise(mean_fdM_san = mean(fdM_san), sd_fdM_san = sd(fdM_san)) %>%
  write_delim("../results/mean_highly_div_windows/mean_fdM_san_diff_fst_f1_parent_windows.txt", delim = "\t")

df_fdM_chr %>% group_by(outlier_q99, outlier_q01) %>%
  dplyr::summarise(mean_fdM = mean(fdM_chr), sd_fdM_chr = sd(fdM_chr)) %>%
  write_delim("../results/mean_highly_div_windows/mean_fdM_chr_diff_fst_f1_parent_windows.txt", delim = "\t")

df_fdM_san_chr %>% group_by(outlier_q99, outlier_q01) %>%
  dplyr::summarise(mean_fdM = mean(fdM_chr), sd_fdM_chr = sd(fdM_chr), mean_fdM_san = mean(fdM_san), sd_fdM_san = sd(fdM_san)) %>%
  write_delim("../results/mean_highly_div_windows/mean_fdM_san_chr_diff_fst_f1_parent_windows.txt", delim = "\t")

### wilcoxon test for difference between group ####

pairwise_pop_interest <- c("CHR_PNG_leuco_BC", "CHR_PNG_leuco_F1", "CHR_PNG_SAN_PNG", "leuco_BC_SAN_PNG", "leuco_F1_SAN_PNG", "leuco_BC_leuco_F1")
pop_interest          <- c("CHR_PNG", "SAN_PNG", "leuco_BC", "leuco_F1")

#### fst ####
df_fst %>%
  filter(pop_Fst %in% pairwise_pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q99") %>%
  group_by(pop_Fst) %>%
  wilcox_test(Fst ~ BG_or_out_diff_F1) %>%
  add_significance() %>%
  write_delim("../results/wilcox_test_highly_div_windows/wilcox_Fst_F1_q01_vs_BG.txt", delim = "\t")

df_fst %>%
  filter(pop_Fst %in% pairwise_pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q99") %>%
  group_by(pop_Fst) %>%
  wilcox_effsize(Fst ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/effsize_wilcox_Fst_F1_q01_vs_BG.txt", delim = "\t")

df_fst %>%
  filter(pop_Fst %in% pairwise_pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q99") %>%
  group_by(pop_Fst) %>%
  cohens_d(Fst ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/cohens_d_wilcox_Fst_F1_q01_vs_BG.txt", delim = "\t")

df_fst %>%
  filter(pop_Fst %in% pairwise_pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q01") %>%
  group_by(pop_Fst) %>%
  rstatix::wilcox_test(Fst ~ BG_or_out_diff_F1) %>%
  add_significance() %>%
  write_delim("../results/wilcox_test_highly_div_windows/wilcox_Fst_F1_q99_vs_BG.txt", delim = "\t")

df_fst %>%
  filter(pop_Fst %in% pairwise_pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q01") %>%
  group_by(pop_Fst) %>%
  wilcox_effsize(Fst ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/effsize_wilcox_Fst_F1_q99_vs_BG.txt", delim = "\t")

df_fst %>%
  filter(pop_Fst %in% pairwise_pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q01") %>%
  group_by(pop_Fst) %>%
  cohens_d(Fst ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/cohens_d_wilcox_Fst_F1_q99_vs_BG.txt", delim = "\t")

#### dxy ####
df_dxy %>%
  filter(pop_dxy %in% pairwise_pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q99") %>%
  group_by(pop_dxy) %>%
  wilcox_test(dxy ~ BG_or_out_diff_F1) %>%
  add_significance() %>%
  write_delim("../results/wilcox_test_highly_div_windows/wilcox_dxy_F1_q01_vs_BG.txt", delim = "\t")

df_dxy %>%
  filter(pop_dxy %in% pairwise_pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q99") %>%
  group_by(pop_dxy) %>%
  wilcox_effsize(dxy ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/effsize_wilcox_dxy_F1_q01_vs_BG.txt", delim = "\t")

df_dxy %>%
  filter(pop_dxy %in% pairwise_pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q99") %>%
  group_by(pop_dxy) %>%
  cohens_d(dxy ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/cohens_d_wilcox_dxy_F1_q01_vs_BG.txt", delim = "\t")

df_dxy %>%
  filter(pop_dxy %in% pairwise_pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q01") %>%
  group_by(pop_dxy) %>%
  wilcox_test(dxy ~ BG_or_out_diff_F1) %>%
  add_significance() %>%
  write_delim("../results/wilcox_test_highly_div_windows/wilcox_dxy_F1_q99_vs_BG.txt", delim = "\t")

df_dxy %>%
  filter(pop_dxy %in% pairwise_pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q01") %>%
  group_by(pop_dxy) %>%
  wilcox_effsize(dxy ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/effsize_wilcox_dxy_F1_q99_vs_BG.txt", delim = "\t")


df_dxy %>%
  filter(pop_dxy %in% pairwise_pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q01") %>%
  group_by(pop_dxy) %>%
  cohens_d(dxy ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/cohens_d_wilcox_dxy_F1_q99_vs_BG.txt", delim = "\t")

#### pi ####
df_pi %>%
  filter(pop_pi %in% pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q99") %>%
  group_by(pop_pi) %>%
  wilcox_test(pi ~ BG_or_out_diff_F1) %>%
  add_significance() %>%
  write_delim("../results/wilcox_test_highly_div_windows/wilcox_pi_F1_q01_vs_BG.txt", delim = "\t")

df_pi %>%
  filter(pop_pi %in% pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q99") %>%
  group_by(pop_pi) %>%
  wilcox_effsize(pi ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/effsize_wilcox_pi_F1_q01_vs_BG.txt", delim = "\t")

df_pi %>%
  filter(pop_pi %in% pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q99") %>%
  group_by(pop_pi) %>%
  cohens_d(pi ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/cohens_d_wilcox_pi_F1_q01_vs_BG.txt", delim = "\t")

df_pi %>%
  filter(pop_pi %in% pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q01") %>%
  group_by(pop_pi) %>%
  wilcox_test(pi ~ BG_or_out_diff_F1) %>%
  add_significance() %>%
  write_delim("../results/wilcox_test_highly_div_windows/wilcox_pi_F1_q99_vs_BG.txt", delim = "\t")

df_pi %>%
  filter(pop_pi %in% pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q01") %>%
  group_by(pop_pi) %>%
  wilcox_effsize(pi ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/effsize_wilcox_pi_F1_q99_vs_BG.txt", delim = "\t")

df_pi %>%
  filter(pop_pi %in% pop_interest) %>%
  filter(!BG_or_out_diff_F1 == "out_q01") %>%
  group_by(pop_pi) %>%
  cohens_d(pi ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/cohens_d_wilcox_pi_F1_q99_vs_BG.txt", delim = "\t")

#### fdM_san #####
df_fdM_san %>%
  filter(!BG_or_out_diff_F1 == "out_q99") %>%
  wilcox_test(fdM_san ~ BG_or_out_diff_F1) %>%
  add_significance() %>%
  write_delim("../results/wilcox_fdM_san_F1_q01_vs_BG.txt", delim = "\t")

df_fdM_san %>%
  filter(!BG_or_out_diff_F1 == "out_q99") %>%
  wilcox_effsize(fdM_san ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/effsize_wilcox_fdM_san_F1_q01_vs_BG.txt", delim = "\t")

df_fdM_san %>%
  filter(!BG_or_out_diff_F1 == "out_q99") %>%
  cohens_d(fdM_san ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/cohens_d_wilcox_fdM_san_F1_q01_vs_BG.txt", delim = "\t")

df_fdM_san %>%
  filter(!BG_or_out_diff_F1 == "out_q01") %>%
  wilcox_test(fdM_san ~ BG_or_out_diff_F1) %>%
  add_significance() %>%
  write_delim("../results/wilcox_test_highly_div_windows/wilcox_fdM_san_F1_q99_vs_BG.txt", delim = "\t")

df_fdM_san %>%
  filter(!BG_or_out_diff_F1 == "out_q01") %>%
  wilcox_effsize(fdM_san ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/effsize_wilcox_fdM_san_F1_q99_vs_BG.txt", delim = "\t")

df_fdM_san %>%
  filter(!BG_or_out_diff_F1 == "out_q01") %>%
  cohens_d(fdM_san ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/cohens_d_wilcox_fdM_san_F1_q99_vs_BG.txt", delim = "\t")

#### fdM_chr #####
df_fdM_chr %>%
  filter(!BG_or_out_diff_F1 == "out_q99") %>%
  wilcox_test(fdM_chr ~ BG_or_out_diff_F1) %>%
  add_significance() %>%
  write_delim("../results/wilcox_fdM_chr_F1_q01_vs_BG.txt", delim = "\t")

df_fdM_chr %>%
  filter(!BG_or_out_diff_F1 == "out_q99") %>%
  wilcox_effsize(fdM_chr ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/effsize_wilcox_fdM_chr_F1_q01_vs_BG.txt", delim = "\t")


df_fdM_chr %>%
  filter(!BG_or_out_diff_F1 == "out_q99") %>%
  cohens_d(fdM_chr ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/cohens_d_wilcox_fdM_chr_F1_q01_vs_BG.txt", delim = "\t")

df_fdM_chr %>%
  filter(!BG_or_out_diff_F1 == "out_q01") %>%
  wilcox_test(fdM_chr ~ BG_or_out_diff_F1) %>%
  add_significance() %>%
  write_delim("../results/wilcox_test_highly_div_windows/wilcox_fdM_chr_F1_q99_vs_BG.txt", delim = "\t")

df_fdM_chr %>%
  filter(!BG_or_out_diff_F1 == "out_q01") %>%
  wilcox_effsize(fdM_chr ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/effsize_wilcox_fdM_chr_F1_q99_vs_BG.txt", delim = "\t")

df_fdM_chr %>%
  filter(!BG_or_out_diff_F1 == "out_q01") %>%
  cohens_d(fdM_chr ~ BG_or_out_diff_F1) %>%
  write_delim("../results/wilcox_test_highly_div_windows/cohens_d_wilcox_fdM_chr_F1_q99_vs_BG.txt", delim = "\t")

#########################################################
#### fdM_san vs fdM_chr in highly divergent windows #####
#########################################################

# put all the fdM values in a single column
df_fdM_long <- pivot_longer(df_fdM_san_chr, cols = c(fdM_chr, fdM_san),
                            names_to = "fdM_test",
                            values_to = "fdM_values")

df_fdM_long %>%
  dplyr::filter(fdM_test %in% c("fdM_chr", "fdM_san")) %>%
  group_by(BG_or_out_diff_F1) %>%
  rstatix::wilcox_test(fdM_values ~ fdM_test, paired = TRUE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  print() %>%
  write_csv2("../results/wilcox_fdm_vs_fdm/wilcox_fdm_chr_vs_fdm_san.csv")

df_fdM_long %>%
  dplyr::filter(fdM_test %in% c("fdM_chr", "fdM_san")) %>%
  group_by(BG_or_out_diff_F1) %>%
  rstatix::cohens_d(fdM_values ~ fdM_test, paired = TRUE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  print() %>%
  write_csv2("../results/wilcox_fdm_vs_fdm/cohens_effect_size_fdm_chr_vs_fdm_san.csv")


### boxplots ####
fdM_san_boxplot_diff_fst_f1_png_q99 <- df_final %>% 
  select(c(scaffold, start, end, outlier_q99, outlier_q01, starts_with("fdM_san"))) %>%
  distinct(scaffold, start, end, outlier_q99, outlier_q01, fdM_san) %>%
  ggplot(aes(x = outlier_q99, y = fdM_san, color = outlier_q99)) +  # Add color mapping here
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.1, alpha = 0.2, width = 0.2) +
  scale_color_manual(values = c("#6C726F", "#41A790")) +  # Ensure color values are set
  labs(y = "fdM_san") +
  coord_cartesian(ylim = c(-0.15, 0.4)) +
  theme(axis.title = element_text(size = 8))

fdM_san_boxplot_diff_fst_f1_png_q01 <- df_final %>%
  select(c(scaffold, start, end, outlier_q99, outlier_q01, starts_with("fdM_san"))) %>%
  distinct(scaffold, start, end, outlier_q99, outlier_q01, fdM_san) %>%
  ggplot(aes(outlier_q01, fdM_san, color = outlier_q01)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size = 0.1, alpha = 0.2, width = 0.2) +
  scale_color_manual(values = c("#6C726F", "#41A790")) +
  labs(y = "fdM_san") +
  coord_cartesian(ylim=c(-0.15, 0.4)) +
  theme(axis.title = element_text(size = 8))

fdM_chr_boxplot_diff_fst_f1_png_q99 <- df_final %>% 
  select(c(scaffold, start, end, outlier_q99, outlier_q01, starts_with("fdM_chr"))) %>%
  distinct(scaffold, start, end, outlier_q99, outlier_q01, fdM_chr) %>%
  ggplot(aes(x = outlier_q99, y = fdM_chr, color = outlier_q99)) +  # Add color mapping here
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.1, alpha = 0.2, width = 0.2) +
  scale_color_manual(values = c("#6C726F", "#8D4F61")) +  # Ensure color values are set
  labs(y = "fdM_chr") +
  coord_cartesian(ylim = c(-0.15, 0.4)) +
  theme(axis.title = element_text(size = 8))

fdM_chr_boxplot_diff_fst_f1_png_q01 <- df_final %>%
  select(c(scaffold, start, end, outlier_q99, outlier_q01, starts_with("fdM_chr"))) %>%
  distinct(scaffold, start, end, outlier_q99, outlier_q01, fdM_chr) %>%
  ggplot(aes(outlier_q01, fdM_chr, color = outlier_q01)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size = 0.1, alpha = 0.2, width = 0.2) +
  scale_color_manual(values = c("#6C726F", "#8D4F61")) +
  labs(y = "fdM_chr") +
  coord_cartesian(ylim=c(-0.15, 0.4)) +
  theme(axis.title = element_text(size = 8))

(fdM_san_boxplot_diff_fst_f1_png_q99 + fdM_san_boxplot_diff_fst_f1_png_q01) / (fdM_chr_boxplot_diff_fst_f1_png_q99 + fdM_chr_boxplot_diff_fst_f1_png_q01)
ggsave("../plots/distribution_fdm_q01_q99.png", units = "cm", width = 30, height = 25)

## fdM, Fst, dxy and pi in backcross windows vs background windows ####

### prepare df ####
df_fst <- df_final %>%
  select(c(scaffold, start, end, bc_sand_window, bc_chryso_window, starts_with("Fst"))) %>%
  pivot_longer(cols = starts_with("Fst"), names_to = "pop_Fst", values_to = "Fst", names_prefix = "Fst_")

df_dxy <- df_final %>%
  select(c(scaffold, start, end, bc_sand_window, bc_chryso_window, starts_with("dxy"))) %>%
  pivot_longer(cols = starts_with("dxy"), names_to = "pop_dxy", values_to = "dxy", names_prefix = "dxy_")

df_pi <- df_final %>%
  select(c(scaffold, start, end, bc_sand_window, bc_chryso_window, starts_with("pi"))) %>%
  pivot_longer(cols = starts_with("pi"), names_to = "pop_pi", values_to = "pi", names_prefix = "pi_")

df_fdM_san <- df_final %>%
  select(c(scaffold, start, end, bc_sand_window, bc_chryso_window, starts_with("fdM_san"))) %>%
  distinct(scaffold, start, end, bc_sand_window, bc_chryso_window,  fdM_san)

df_fdM_chr <- df_final %>%
  select(c(scaffold, start, end, bc_sand_window, bc_chryso_window, starts_with("fdM_chr"))) %>%
  distinct(scaffold, start, end, bc_sand_window, bc_chryso_window,  fdM_chr)

### extract mean values ####
df_fst %>% group_by(bc_sand_window, pop_Fst) %>%
  dplyr::summarise(mean_Fst = mean(Fst), sd_Fst = sd(Fst)) %>%
  write_delim("../results/mean_backcross_windows/mean_Fst_bc_sand_windows.txt", delim = "\t")

df_dxy %>% group_by(bc_sand_window, pop_dxy) %>%
  dplyr::summarise(mean_dxy = mean(dxy), sd_dxy = sd(dxy)) %>%
  write_delim("../results/mean_backcross_windows/mean_dxy_bc_sand_windows.txt", delim = "\t")

df_pi %>% group_by(bc_sand_window, pop_pi) %>%
  dplyr::summarise(mean_pi = mean(pi), sd_pi = sd(pi)) %>%
  write_delim("../results/mean_backcross_windows/mean_pi_bc_sand_windows.txt", delim = "\t")

df_fdM_san %>% group_by(bc_sand_window) %>%
  dplyr::summarise(mean_fdM_san = mean(fdM_san), sd_fdM_san = sd(fdM_san)) %>%
  write_delim("../results/mean_backcross_windows/mean_fdM_san_bc_sand_windows.txt", delim = "\t")

df_fdM_chr %>% group_by(bc_sand_window) %>%
  dplyr::summarise(mean_fdM_chr = mean(fdM_chr), sd_fdM_chr = sd(fdM_chr)) %>%
  write_delim("../results/mean_backcross_windows/mean_fdM_chr_bc_sand_windows.txt", delim = "\t")


### wilcoxon test for difference between groups ####

pairwise_pop_interest <- c("CHR_PNG_leuco_BC", "CHR_PNG_leuco_F1", "CHR_PNG_SAN_PNG", "leuco_BC_SAN_PNG", "leuco_F1_SAN_PNG", "leuco_BC_leuco_F1")
pop_interest          <- c("CHR_PNG", "SAN_PNG", "leuco_BC", "leuco_F1")

#### fst ####
df_fst %>%
  filter(pop_Fst %in% pairwise_pop_interest) %>%
  group_by(pop_Fst) %>%
  wilcox_test(Fst ~ bc_sand_window) %>%
  add_significance() %>%
  write_delim("../results/wilcox_test_backcross_windows/wilcox_fst_BC_vs_BG.txt", delim = "\t")

df_fst %>%
  filter(pop_Fst %in% pairwise_pop_interest) %>%
  group_by(pop_Fst) %>%
  wilcox_effsize(Fst ~ bc_sand_window) %>%
  write_delim("../results/wilcox_test_backcross_windows/effsize_wilcox_fst_BC_vs_BG.txt", delim = "\t")

df_fst %>%
  filter(pop_Fst %in% pairwise_pop_interest) %>%
  group_by(pop_Fst) %>%
  cohens_d(Fst ~ bc_sand_window) %>%
  write_delim("../results/wilcox_test_backcross_windows/cohens_d_wilcox_fst_BC_vs_BG.txt", delim = "\t")

#### dxy ####
df_dxy %>%
  filter(pop_dxy %in% pairwise_pop_interest) %>%
  group_by(pop_dxy) %>%
  wilcox_test(dxy ~ bc_sand_window) %>%
  add_significance() %>%
  write_delim("../results/wilcox_test_backcross_windows/wilcox_dxy_BC_vs_BG.txt", delim = "\t")

df_dxy %>%
  filter(pop_dxy %in% pairwise_pop_interest) %>%
  group_by(pop_dxy) %>%
  wilcox_effsize(dxy ~ bc_sand_window) %>%
  write_delim("../results/wilcox_test_backcross_windows/effsize_wilcox_dxy_BC_vs_BG.txt", delim = "\t")

df_dxy %>%
  filter(pop_dxy %in% pairwise_pop_interest) %>%
  group_by(pop_dxy) %>%
  cohens_d(dxy ~ bc_sand_window) %>%
  write_delim("../results/wilcox_test_backcross_windows/cohens_d_wilcox_dxy_BC_vs_BG.txt", delim = "\t")


#### pi ####
df_pi %>%
  filter(pop_pi %in% pop_interest) %>%
  group_by(pop_pi) %>%
  wilcox_test(pi ~ bc_sand_window) %>%
  add_significance() %>%
  write_delim("../results/wilcox_test_backcross_windows/wilcox_pi_BC_vs_BG.txt", delim = "\t")

df_pi %>%
  filter(pop_pi %in% pop_interest) %>%
  group_by(pop_pi) %>%
  wilcox_effsize(pi ~ bc_sand_window) %>%
  write_delim("../results/wilcox_test_backcross_windows/effsize_wilcox_pi_BC_vs_BG.txt", delim = "\t")

df_pi %>%
  filter(pop_pi %in% pop_interest) %>%
  group_by(pop_pi) %>%
  cohens_d(pi ~ bc_sand_window) %>%
  write_delim("../results/wilcox_test_backcross_windows/cohens_d_wilcox_pi_BC_vs_BG.txt", delim = "\t")

#### fdM_san #####
df_fdM_san %>%
  wilcox_test(fdM_san ~ bc_sand_window) %>%
  add_significance() %>%
  write_delim("../results/wilcox_test_backcross_windows/wilcox_fdM_san_BC_vs_BG.txt", delim = "\t")

df_fdM_san %>%
  wilcox_effsize(fdM_san ~ bc_sand_window) %>%
  write_delim("../results/wilcox_test_backcross_windows/effsize_wilcox_fdM_san_BC_vs_BG.txt", delim = "\t")

df_fdM_san %>%
  cohens_d(fdM_san ~ bc_sand_window) %>%
  write_delim("../results/wilcox_test_backcross_windows/cohens_d_wilcox_fdM_san_BC_vs_BG.txt", delim = "\t")

#### fdM_chr #####
df_fdM_chr %>%
  wilcox_test(fdM_chr ~ bc_sand_window) %>%
  add_significance() %>%
  write_delim("../results/wilcox_test_backcross_windows/wilcox_fdM_chr_BC_vs_BG.txt", delim = "\t")

df_fdM_chr %>%
  wilcox_effsize(fdM_chr ~ bc_sand_window) %>%
  write_delim("../results/wilcox_test_backcross_windows/effsize_wilcox_fdM_chr_BC_vs_BG.txt", delim = "\t")

df_fdM_chr %>%
  cohens_d(fdM_chr ~ bc_sand_window) %>%
  write_delim("../results/wilcox_test_backcross_windows/cohens_d_wilcox_fdM_chr_BC_vs_BG.txt", delim = "\t")


## barrier loci among parental species ####

### for bc sand windows ####
table_barrier_bc <- df_final %>%
  select(c(scaffold, start, end, bc_sand_window, bc_chryso_window, barrier_loci)) %>%
  distinct(scaffold, start, end, bc_sand_window, bc_chryso_window, barrier_loci)  %>%
  tabyl(bc_sand_window, barrier_loci) #create directly a table we can use for fisher test

janitor::chisq.test(table_barrier_bc)

table_barrier_bc %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() %>%
  write_delim("../results/chisq_barrier_loci_bc_sand.txt", delim = "\t")

### for diff fst F1-parents ####
table_barrier_fst_f1_q99 <- df_final %>%
  select(c(scaffold, start, end, outlier_q99, outlier_q01, barrier_loci)) %>%
  distinct(scaffold, start, end, outlier_q99, outlier_q01, barrier_loci)  %>%
  tabyl(outlier_q99, barrier_loci) #create directly a table we can use for fisher test

table_barrier_fst_f1_q01 <- df_final %>%
  select(c(scaffold, start, end, outlier_q99, outlier_q01, barrier_loci)) %>%
  distinct(scaffold, start, end, outlier_q99, outlier_q01, barrier_loci)  %>%
  tabyl(outlier_q01, barrier_loci) #create directly a table we can use for fisher test

janitor::fisher.test(table_barrier_fst_f1_q99)
janitor::chisq.test(table_barrier_fst_f1_q01)

table_barrier_fst_f1_q01 %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() %>%
  write_delim("../results/chisq_barrier_loci_diff_fst_f1_parent_q01.txt", delim = "\t")

table_barrier_fst_f1_q99 %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() %>%
  write_delim("../results/chisq_barrier_loci_diff_fst_f1_parent_q99.txt", delim = "\t")



#############################################################################|
# Are the outlier windows randomely distributed across chromosome or not? ####
#############################################################################|

## Number of backcross windows per chromosome ####

df_final$scaffold    = factor(df_final$scaffold)
df_final$pop         = factor(df_final$pop)
df_final$bc_sand_window = factor(df_final$bc_sand_window)
nb_bc_wind <- df_final %>%
  group_by(scaffold, bc_sand_window, .drop=FALSE) %>%
  tally %>% 
  dplyr::filter(bc_sand_window == "bc_sand_win") %>%
  dplyr::select(c(scaffold, nb_bc_win = n)) %>%
  ungroup()

## Number of outlier diff-fst windows for F1-parents comparison ####

df_final$scaffold    = factor(df_final$scaffold)
df_final$pop         = factor(df_final$pop)
df_final$outlier_q99 = factor(df_final$outlier_q99)
nb_outlier_win_chr_fst <- df_final %>%
  group_by(pop, scaffold, outlier_q99, .drop=FALSE) %>%
  tally %>% 
  dplyr::filter(outlier_q99 == "outlier") %>%
  dplyr::select(c(pop, scaffold, nb_outlier_win = n)) %>%
  ungroup()
  

## Reshuffling windows ####

### function to reshuffle a column from a dataset ####

shuffle <- function(.data, n, perm_cols){
  cols_ids <- match(perm_cols, colnames(.data))
  #creates a sequence of integers from 1 to the number of rows of the dataset
  ids <- seq_len(nrow(.data))
  #generates n permutations of the vector of IDs
  n_ids <- rerun(n, sample(ids))
  
  map_dfr(n_ids, function(x){
    .data[ids, cols_ids] <- .data[x, cols_ids]
    .data
  })
}

### function to do the permutation and calculate p-value for diff-fst outlier windows ####

fun_perm_pvalue <- function(diff_dataset, nb_outlier_dataset, population, statistics) {
  set.seed(123)
  perm_df <- diff_dataset %>%
    filter(pop == population) %>%
    dplyr::select(c(scaffold, start, end, pop, outlier_q99)) %>%
    shuffle(n = 10000, perm_cols = c("outlier_q99")) %>%
    mutate(perm = rep(1:10000, length.out = n()))  
  nb_outlier_dataset <- nb_outlier_dataset %>%
    mutate(scaffold = as.numeric(scaffold))
  
  
  #### extract information from the permutation (mean across chromosome and permutation)
  perm_df_stat <- perm_df %>%
    group_by(scaffold, perm, outlier_q99, .drop=FALSE) %>%
    tally %>%
    dplyr::filter(outlier_q99 == "outlier") %>%
    dplyr::select(c(perm, scaffold, nb_outlier_win_exp = n)) %>%
    ungroup() %>%
    group_by(scaffold) %>%
    dplyr::mutate(mean_nb_outlier_win_exp = mean(nb_outlier_win_exp)) %>%
    dplyr::mutate(sd_nb_outlier_win_exp = sd(nb_outlier_win_exp)) %>%
    ungroup() %>%
    left_join(select(filter(nb_outlier_dataset, pop == population), c(scaffold, nb_outlier_win_obs = nb_outlier_win)), c("scaffold" = "scaffold")) %>%
    mutate(scaffold = as.character(scaffold))
  
  #### calculate pvalue
  pvalue_df <- perm_df_stat %>%
    group_by(scaffold) %>%
    mutate(pvalue = sum(abs(nb_outlier_win_exp) >= abs(nb_outlier_win_obs)) / 10000) %>% #change according to number of permutation
    filter(perm == 1) %>%
    select(scaffold, mean_nb_outlier_win_exp, sd_nb_outlier_win_exp, nb_outlier_win_obs, pvalue) %>%
    mutate(pvalue_corr = p.adjust(pvalue, method = "BH", n = 24))
  
  #### save pvalue table
  write.table(pvalue_df, file = paste("../results/permutation_test/pvalue_out_window_", population, "_", statistics, ".txt", sep = ""), sep = "\t", 
              row.names = FALSE)
  
}

### apply function to each population and statistics ####

#### outlier diff-fst windows between F1 and parents #####
fun_perm_pvalue(df_final, nb_outlier_win_chr_fst, "f1_png", "fst")
fun_perm_pvalue(df_final, nb_outlier_win_chr_fst, "f1_slm", "fst")
fun_perm_pvalue(df_final, nb_outlier_win_chr_fst, "f1_allo", "fst")


### function to do the permutation and calculate p-value for BC-sand outlier windows ####

fun_perm_BC_pvalue <- function(dataset, nb_outlier_dataset) {
  set.seed(123)
  perm_df <- dataset %>%
    dplyr::select(c(scaffold, start, end, bc_sand_window)) %>%
    group_by(scaffold) %>%
    distinct(start, end, bc_sand_window) %>%
    ungroup() %>%
    shuffle(n = 10000, perm_cols = c("bc_sand_window")) %>%
    mutate(perm = rep(c(1:10000), length.out = n()))
  
  nb_outlier_dataset <- nb_outlier_dataset %>%
    mutate(scaffold = as.numeric(scaffold))
  
  #### extract information from the permutation (mean across chromosome and permutation)
  perm_df_stat <- perm_df %>%
    group_by(scaffold, perm, bc_sand_window, .drop=FALSE) %>%
    tally %>%
    dplyr::filter(bc_sand_window == "bc_sand_win") %>%
    dplyr::select(c(perm, scaffold, nb_outlier_win_exp = n)) %>%
    ungroup() %>%
    group_by(scaffold) %>%
    dplyr::mutate(mean_nb_outlier_win_exp = mean(nb_outlier_win_exp)) %>%
    dplyr::mutate(sd_nb_outlier_win_exp = sd(nb_outlier_win_exp)) %>%
    ungroup() %>%
    left_join(select(nb_outlier_dataset, c(scaffold, nb_outlier_win_obs = nb_bc_win)), by = "scaffold") #change nb_bc_win according to name of the column in the outlier windows input
  
  #### calculate pvalue
  pvalue_df <- perm_df_stat %>%
    group_by(scaffold) %>%
    mutate(pvalue = sum(abs(nb_outlier_win_exp) >= abs(nb_outlier_win_obs)) / 10000) %>% #change according to number of permutation
    filter(perm == 1) %>%
    select(scaffold, mean_nb_outlier_win_exp, sd_nb_outlier_win_exp, nb_outlier_win_obs, pvalue) %>%
    mutate(pvalue_corr = p.adjust(pvalue, method = "BH", n = 24))
  
  #### save pvalue table
  write.table(pvalue_df, file = "../results/permutation_test/pvalue_bc_sand_window.txt",  sep = "\t", row.names = FALSE)
  
}

#### bc-sand windows ####
fun_perm_BC_pvalue(df_final, nb_bc_wind)

  
#############################################################################|
# Is the difference in fst and dxy significantly different from zero?    ####
#############################################################################|

## function to compute wilcoxon-test across chromosome and for the whole-genome (scaffold = 25) ####

wilcox_diff <- function(dataset, pop, statistic){
  df1 <- dataset %>%
    dplyr::filter(pop == pop) %>%
    dplyr::select(scaffold, start, diff = paste(statistic, "_diff", sep = ""))
  pvalue_list  <- list()
  #wilcoxon test for each chromosome
  for (i in 1:24){
    wilcox_df <- tibble(scaffold = double(), pvalue = double(), v_stat = double())
    df2       <- dplyr::filter(df1, scaffold == i)
    wilcox    <- wilcox.test(df2$diff, mu = 0)
    p <- wilcox$p.value
    v <- as.vector(wilcox$statistic)
    wilcox_df <- tibble(scaffold = i, pvalue = p, v_stat = v)
    pvalue_list[[i]] <- wilcox_df
  }
  #WG wilcoxon test 
  wilcox_df <- tibble(scaffold = double(), pvalue = double(), v_stat = double())
  wg_wilcox <- wilcox.test(df1$diff, mu = 0)
  p <- wg_wilcox$p.value
  v <- as.vector(wg_wilcox$statistic)
  wilcox_df <- tibble(scaffold = 25, pvalue = p, v_stat = v)
  pvalue_list[[25]] <- wilcox_df
  pvalue_wilcox <- dplyr::bind_rows(pvalue_list)
  write_delim(pvalue_wilcox, file = paste("../results/wilcox_diff_from_zero/wilcox_", pop, "_", statistic, ".txt", sep = ""), delim = " ")
}

## run the function for each population and statistics ####
### /!\ when running the script with the function, we have always the same result. But if
### we run the script outside the function it's working. So that's what we did in the end

wilcox_diff(df_final, "f1_allo", "dxy")
wilcox_diff(df_final, "f1_png", "dxy")
wilcox_diff(df_final, "f1_slm", "dxy")
wilcox_diff(df_final, "bc_allo", "dxy")
wilcox_diff(df_final, "bc_png", "dxy")
wilcox_diff(df_final, "bc_slm", "dxy")

wilcox_diff(df_final, "f1_allo", "fst")
wilcox_diff(df_final, "f1_png", "fst")
wilcox_diff(df_final, "f1_slm", "fst")
wilcox_diff(df_final, "bc_allo", "fst")
wilcox_diff(df_final, "bc_png", "fst")
wilcox_diff(df_final, "bc_slm", "fst")


### is there a bias towards some chromosme for the backcross windows ####
df_bc <- df_final %>%
  filter(pop == "bc_png") %>%
  dplyr::select(scaffold, start, bc_sand_window, pop) %>%
  #filter(bc_sand_window == "bc_sand_win") %>%
  group_by(scaffold, bc_sand_window, .drop = FALSE) %>%
  #distinct(scaffold, start, bc_sand_window) %>%
  tally() %>%
  filter(bc_sand_window == "bc_sand_win") %>%
  select(!bc_sand_window) %>%
  dplyr::rename("nb_bc_wind" = n) %>%
  ungroup() %>%
  #add manually the chromosomes with 0 windows which were removed
  add_row(scaffold = c(3, 6, 7, 8, 10, 13, 15, 18), nb_bc_wind = rep(0, 8))

krt <- kruskal.test(nb_bc_wind~as.factor(scaffold), data = df_bc)
pva <- krt$p.value
chi <- as.vector(krt$statistic)
pwt <- pairwise.wilcox.test(df_bc$nb_bc_wind, as.factor(df_bc$scaffold), p.adjust.method = "BH")
write.table(pwt$p.value, file = "../results/wilcox_test_backcross_windows/pvalue_nb_bc_wind_per_chromosome.txt", row.names = TRUE)



