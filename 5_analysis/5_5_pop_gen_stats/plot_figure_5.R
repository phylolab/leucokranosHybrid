###############################|
## 
## Plotting figure 5
## A. leucokranos
##
## author: Sarah Schmid
##
## june 2022
################################|


## Load libraries ####
library(ggrepel)
library(tidyverse)
library(patchwork)
library(here)

# plot Fst, dxy, fdM ####

col.barrier <- c(barrier = "#FFA17A", no = "#7e7f80")
col.density <- c("#BDC581", "#25CCF7", "#FC427B")
col.density.2 <- c("Allopatric" = "#BDC581", "Papua New Guinea" = "#25CCF7", "Solomon Islands" = "#FC427B")
theme_set(theme_light(base_size = 12))

## load data ####

### barrier loci
barrier_loci    <- read_delim("../results/outlier_fst_parents.txt") %>%
  mutate(barrier_loci = rep("barrier")) %>%
  select(c(scaffold, start, barrier_loci))

### fdm 
data_fdm_san <- read_csv("../results/fd.san_aus.san_png.chr_png.cla.w50m100.csv.gz") %>% 
  select(scaffold, start, end, fdM) %>%
  rename(fdM_san = fdM)
data_fdm_chr <- read_csv("../results/fd.chr_fiji.chr_png.san_png.cla.w50m100.csv.gz") %>% 
  select(scaffold, start, end, fdM) %>%
  rename(fdM_chr = fdM)

### combine everything
fig5_df <- read_csv("../results/all_species_popgen.w50m100.csv.gz") %>%
  dplyr::select(c(scaffold, start, end, pi_SAN_PNG, pi_SAN_AUS, pi_SAN_SLM, pi_CHR_PNG, pi_CHR_FIJI, pi_CHR_SLM, ends_with(c("CHR_SLM_SAN_SLM", "CHR_FIJI_SAN_AUS", "CHR_PNG_SAN_PNG")))) %>%
  left_join(data_fdm_san) %>%
  left_join(data_fdm_chr) %>%
  left_join(barrier_loci) %>%
  mutate(barrier_loci = replace_na(barrier_loci, "no")) %>%
  mutate(scaffold = as.numeric(str_remove(scaffold, "^Chr")))


## plotting dxy, fst, pi, fdm ####

### fst ####
mean_fst <- fig5_df %>%
  dplyr::filter(!is.na(Fst_CHR_PNG_SAN_PNG)) %>%  # Remove rows where Fst_CHR_PNG_SAN_PNG is NA
  select(c(scaffold, start, end, Fst_CHR_PNG_SAN_PNG, barrier_loci)) %>%
  group_by(scaffold) %>%
  dplyr::summarise(mean_fst = mean(Fst_CHR_PNG_SAN_PNG)) %>%
  ungroup()

df_fst <- fig5_df %>%
  select(c(scaffold, start, end, Fst_CHR_PNG_SAN_PNG, barrier_loci)) %>%
  pivot_longer(cols = starts_with("Fst"), names_to = "pop_fst", values_to = "Fst", names_prefix = "Fst_") %>%
  mutate(q_99 = quantile(Fst, probs = 0.99, na.rm = T)) %>% 
  mutate(q_01 = quantile(Fst, probs = 0.01, na.rm = T)) %>%
  mutate(outlier_q99  = ifelse(Fst >= q_99 | Fst <= q_01, "outlier", "no"))

plot_fst <- ggplot(df_fst, aes(start, Fst)) + 
  geom_rect(data = dplyr::filter(df_fst, scaffold %% 2 == 1),
            fill = "grey92", color = "grey92", 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + #color only odd chromosomes panel
  geom_point(color = "#7e7f80") +
  geom_point(data = dplyr::filter(df_fst, barrier_loci == "barrier"), aes(start, Fst), color = "#FFA17A") +
  #scale_color_manual(values = col.barrier) +
  #geom_hline(data = diff_df, aes(yintercept = q_99), size = 0.5, linetype= "dotted",  color = "black") +
  #geom_hline(data = diff_df, aes(yintercept = q_01), size = 0.5, linetype= "dotted",  color = "black") +
  geom_hline(data = mean_fst, aes(yintercept = mean_fst), linewidth = 0.4, color = "#464747", linetype = "dotted") +
  facet_grid(~scaffold, scales = 'free_x', space = 'free_x', switch = 'x', labeller = labeller(scaffold = c(1:24))) +
  labs(y = expression(F[ST])) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank(),
        #strip.text.x = element_text(color = "black"),
        strip.text.x = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.background = element_blank())

### dxy ####
mean_dxy <- fig5_df %>% 
  dplyr::filter(!is.na(dxy_CHR_PNG_SAN_PNG)) %>%  # Remove rows where dxy_CHR_PNG_SAN_PNG is NA
  select(c(scaffold, start, end, dxy_CHR_PNG_SAN_PNG, barrier_loci)) %>% 
  group_by(scaffold) %>% 
  dplyr::summarise(mean_dxy = mean(dxy_CHR_PNG_SAN_PNG)) %>% 
  ungroup()

df_dxy <- fig5_df %>%
  select(c(scaffold, start, end, dxy_CHR_PNG_SAN_PNG, barrier_loci)) %>%
  pivot_longer(cols = starts_with("dxy"), names_to = "pop_dxy", values_to = "dxy", names_prefix = "dxy_") %>%
  mutate(q_99 = quantile(dxy, probs = 0.99, na.rm = T)) %>% 
  mutate(q_01 = quantile(dxy, probs = 0.01, na.rm = T)) %>%
  mutate(outlier_q99  = ifelse(dxy >= q_99 | dxy <= q_01, "outlier", "no"))

plot_dxy <- ggplot(df_dxy, aes(start, dxy)) + 
  geom_rect(data = dplyr::filter(df_dxy, scaffold %% 2 == 1),
            fill = "grey92", color = "grey92", 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + #color only odd chromosomes panel
  geom_point(color = "#7e7f80") +
  geom_point(data = dplyr::filter(df_dxy, barrier_loci == "barrier"), aes(start, dxy), color = "#FFA17A") +
  #scale_color_manual(values = col.barrier) +
  #geom_hline(data = diff_df, aes(yintercept = q_99), size = 0.5, linetype= "dotted",  color = "black") +
  #geom_hline(data = diff_df, aes(yintercept = q_01), size = 0.5, linetype= "dotted",  color = "black") +
  geom_hline(data = mean_dxy, aes(yintercept = mean_dxy), linewidth = 0.4, color = "#464747", linetype = "dotted") +
  facet_grid(~scaffold, scales = 'free_x', space = 'free_x', switch = 'x', labeller = labeller(scaffold = c(1:24))) +
  labs(y = expression(d[XY])) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank(),
        #strip.text.x = element_text(color = "black"),
        strip.text.x = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.background = element_blank())

### pi san ####
mean_pi_san <- fig5_df %>%
  dplyr::filter(!is.na(pi_SAN_PNG)) %>%  # Remove rows where pi_SAN_PNG is NA
  select(c(scaffold, start, end, pi_SAN_PNG, barrier_loci)) %>%
  group_by(scaffold) %>%
  dplyr::summarise(mean_pi_san = mean(pi_SAN_PNG, na.rm = TRUE)) %>%
  ungroup()

df_pi_san <- fig5_df %>%
  select(c(scaffold, start, end, pi_SAN_PNG, barrier_loci)) %>%
  pivot_longer(cols = starts_with("pi"), names_to = "pop_pi", values_to = "pi", names_prefix = "pi_") %>%
  mutate(q_99 = quantile(pi, probs = 0.99, na.rm = T)) %>% 
  mutate(q_01 = quantile(pi, probs = 0.01, na.rm = T)) %>%
  mutate(outlier_q99  = ifelse(pi >= q_99 | pi <= q_01, "outlier", "no"))

plot_pi_san <- ggplot(df_pi_san, aes(start, pi)) + 
  geom_rect(data = dplyr::filter(df_pi_san, scaffold %% 2 == 1),
            fill = "grey92", color = "grey92", 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + #color only odd chromosomes panel
  geom_point(color = "#7e7f80") +
  geom_point(data = dplyr::filter(df_pi_san, barrier_loci == "barrier"), aes(start, pi), color = "#FFA17A") +
  #scale_color_manual(values = col.barrier) +
  #geom_hline(data = diff_df, aes(yintercept = q_99), size = 0.5, linetype= "dotted",  color = "black") +
  #geom_hline(data = diff_df, aes(yintercept = q_01), size = 0.5, linetype= "dotted",  color = "black") +
  geom_hline(data = mean_pi_san, aes(yintercept = mean_pi_san), linewidth = 0.4, color = "#464747", linetype = "dotted") +
  facet_grid(~scaffold, scales = 'free_x', space = 'free_x', switch = 'x', labeller = labeller(scaffold = c(1:24))) +
  labs(y = expression(italic("π A. sandaracinos"))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank(),
        #strip.text.x = element_text(color = "black"),
        strip.text.x = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.background = element_blank())

### pi chryso ####
mean_pi_chr <- fig5_df %>%
  dplyr::filter(!is.na(pi_CHR_PNG)) %>%  # Remove rows where pi_CHR_PNG is NA
  select(c(scaffold, start, end, pi_CHR_PNG, barrier_loci)) %>%
  group_by(scaffold) %>%
  dplyr::summarise(mean_pi_chr = mean(pi_CHR_PNG)) %>%
  ungroup()

df_pi_chr <- fig5_df %>%
  select(c(scaffold, start, end, pi_CHR_PNG, barrier_loci)) %>%
  pivot_longer(cols = starts_with("pi"), names_to = "pop_pi", values_to = "pi", names_prefix = "pi_") %>%
  mutate(q_99 = quantile(pi, probs = 0.99, na.rm = T)) %>% 
  mutate(q_01 = quantile(pi, probs = 0.01, na.rm = T)) %>%
  mutate(outlier_q99  = ifelse(pi >= q_99 | pi <= q_01, "outlier", "no"))

plot_pi_chr <- ggplot(df_pi_chr, aes(start, pi)) + 
  geom_rect(data = dplyr::filter(df_pi_chr, scaffold %% 2 == 1),
            fill = "grey92", color = "grey92", 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + #color only odd chromosomes panel
  geom_point(color = "#7e7f80") +
  geom_point(data = dplyr::filter(df_pi_chr, barrier_loci == "barrier"), aes(start, pi), color = "#FFA17A") +
  #scale_color_manual(values = col.barrier) +
  #geom_hline(data = diff_df, aes(yintercept = q_99), size = 0.5, linetype= "dotted",  color = "black") +
  #geom_hline(data = diff_df, aes(yintercept = q_01), size = 0.5, linetype= "dotted",  color = "black") +
  geom_hline(data = mean_pi_chr, aes(yintercept = mean_pi_chr), linewidth = 0.4, color = "#464747", linetype = "dotted") +
  facet_grid(~scaffold, scales = 'free_x', space = 'free_x', switch = 'x', labeller = labeller(scaffold = c(1:24))) +
  labs(y = expression(italic("pi A. chrysopterus"))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank(),
        #strip.text.x = element_text(color = "black"),
        strip.text.x = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.background = element_blank())

### fdm san ####
mean_fdM_san <- fig5_df %>%
  drop_na() %>%
  select(c(scaffold, start, end, fdM_san, barrier_loci)) %>%
  group_by(scaffold) %>%
  dplyr::summarise(mean_fdm_san = mean(fdM_san, na.rm = TRUE)) %>%
  ungroup()

df_fdM_san <- fig5_df %>% 
  select(c(scaffold, start, end, fdM_san, barrier_loci)) %>% 
  pivot_longer(cols = starts_with("fdM_san"), names_to = "pop_fdM_san", values_to = "fdM_san", names_prefix = "fdM_san") %>% 
  mutate(q_99 = quantile(fdM_san, probs = 0.99, na.rm = T)) %>% 
  mutate(q_01 = quantile(fdM_san, probs = 0.01, na.rm = T)) %>% 
  mutate(outlier_q99 = ifelse(fdM_san >= q_99 | fdM_san <= q_01, "outlier", "no")) %>%
  left_join(mean_fdM_san, by = "scaffold")


plot_fdM_san <- ggplot(df_fdM_san, aes(start, fdM_san)) +  
  geom_rect(data = dplyr::filter(df_fdM_san, scaffold %% 2 == 1),
            fill = "grey92", color = "grey92",  
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + # Color only odd chromosomes panel
  geom_point(color = "#7e7f80") +
  geom_point(data = dplyr::filter(df_fdM_san, barrier_loci == "barrier"), aes(start, fdM_san), color = "#FFA17A") +
  geom_hline(aes(yintercept = mean_fdm_san), linewidth = 0.4, color = "#464747", linetype = "dotted") +  # Using the mean from df_fdM_san
  facet_grid(~scaffold, scales = 'free_x', space = 'free_x', switch = 'x', labeller = labeller(scaffold = c(1:24))) +
  labs(y = expression(italic(F)[dM]("A. sand"))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank(),
        #strip.text.x = element_text(color = "black"),
        strip.text.x = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.background = element_blank())


### fdm chr ####
mean_fdM_chr <- fig5_df %>%
  drop_na() %>%
  select(c(scaffold, start, end, fdM_chr, barrier_loci)) %>%
  group_by(scaffold) %>%
  dplyr::summarise(mean_fdm_chr = mean(fdM_chr, na.rm = TRUE)) %>%
  ungroup()

df_fdM_chr <- fig5_df %>% 
  select(c(scaffold, start, end, fdM_chr, barrier_loci)) %>% 
  pivot_longer(cols = starts_with("fdM_chr"), names_to = "pop_fdM_chr", values_to = "fdM_chr", names_prefix = "fdM_chr") %>% 
  mutate(q_99 = quantile(fdM_chr, probs = 0.99, na.rm = T)) %>% 
  mutate(q_01 = quantile(fdM_chr, probs = 0.01, na.rm = T)) %>% 
  mutate(outlier_q99 = ifelse(fdM_chr >= q_99 | fdM_chr <= q_01, "outlier", "no")) %>%
  left_join(mean_fdM_chr, by = "scaffold")


plot_fdM_chr <- ggplot(df_fdM_chr, aes(start, fdM_chr)) +  
  geom_rect(data = dplyr::filter(df_fdM_chr, scaffold %% 2 == 1),
            fill = "grey92", color = "grey92",  
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + # Color only odd chromosomes panel
  geom_point(color = "#7e7f80") +
  geom_point(data = dplyr::filter(df_fdM_chr, barrier_loci == "barrier"), aes(start, fdM_chr), color = "#FFA17A") +
  geom_hline(aes(yintercept = mean_fdm_chr), linewidth = 0.4, color = "#464747", linetype = "dotted") +  # Using the mean from df_fdM_chr
  facet_grid(~scaffold, scales = 'free_x', space = 'free_x', switch = 'x', labeller = labeller(scaffold = c(1:24))) +
  labs(y = expression(italic(F)[dM]("A. chryso"))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank(),
        strip.text.x = element_text(color = "black"),
        panel.border = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title = element_text(size = 8),
        axis.l = element_blank(),
        panel.background = element_blank())

## plotting density among populations ####

### fst ####

fst_df <-  fig5_df %>%
  select(c(scaffold, start, end, Fst_CHR_FIJI_SAN_AUS, Fst_CHR_PNG_SAN_PNG, Fst_CHR_SLM_SAN_SLM)) %>%
  dplyr::rename(Fst_Allopatric = Fst_CHR_FIJI_SAN_AUS) %>%
  dplyr::rename(Fst_Papua_New_Guinea = Fst_CHR_PNG_SAN_PNG) %>%
  dplyr::rename(Fst_Solomon_Islands = Fst_CHR_SLM_SAN_SLM) %>%
  pivot_longer(cols = starts_with("Fst_"), names_to = "pairwise_comb", values_to = "fst_value", names_prefix = "Fst_")

fst_labels <- fst_df %>% 
  group_by(pairwise_comb) %>%
  dplyr::summarise(
    xpos = max(density(fst_value, na.rm = TRUE)$x[which.max(density(fst_value, na.rm = TRUE)$y)]),
    ypos = max(density(fst_value, na.rm = TRUE)$y))

fst_density_WG <- fst_df %>%
  ggplot(aes(fst_value, color = pairwise_comb)) + 
  geom_density(bw = 0.02) +
  coord_cartesian(xlim = c(0.25, 1)) +
  scale_color_manual(values = col.density) +
  labs(x = expression(F[ST]), y = "Density") +
  geom_text_repel(data = fst_labels, aes(x = xpos, y = ypos, label = pairwise_comb), nudge_y = 1) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8))


### dxy ####

dxy_df <-  fig5_df %>%
  select(c(scaffold, start, end, dxy_CHR_FIJI_SAN_AUS, dxy_CHR_PNG_SAN_PNG, dxy_CHR_SLM_SAN_SLM)) %>%
  dplyr::rename(dxy_Allopatric = dxy_CHR_FIJI_SAN_AUS) %>%
  dplyr::rename(dxy_Papua_New_Guinea = dxy_CHR_PNG_SAN_PNG) %>%
  dplyr::rename(dxy_Solomon_Islands = dxy_CHR_SLM_SAN_SLM) %>%
  pivot_longer(cols = starts_with("dxy_"), names_to = "pairwise_comb", values_to = "dxy_value", names_prefix = "dxy_")

dxy_labels <- dxy_df %>% 
  group_by(pairwise_comb) %>%
  dplyr::summarise(
    xpos = max(density(dxy_value, na.rm = TRUE)$x[which.max(density(dxy_value, na.rm = TRUE)$y)]),
    ypos = max(density(dxy_value, na.rm = TRUE)$y))

dxy_density_WG <- dxy_df %>%
  ggplot(aes(dxy_value, color = pairwise_comb)) + 
  geom_density(bw = 0.02) +
  coord_cartesian(xlim = c(0.25, 1)) +
  scale_color_manual(values = col.density) +
  labs(x = expression(d[XY]), y = "Density") +
  geom_text_repel(data = dxy_labels, aes(x = xpos, y = ypos, label = pairwise_comb), nudge_y = 1) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8))


### pi ####

pi_df <-  fig5_df %>%
  select(c(scaffold, start, end, pi_CHR_FIJI, pi_SAN_AUS, pi_CHR_PNG, pi_SAN_PNG, pi_CHR_SLM, pi_SAN_SLM)) %>%
  dplyr::rename(pi_CHP_Fiji = pi_CHR_FIJI) %>%
  dplyr::rename(pi_SAN_Christmas_Ids = pi_SAN_AUS) %>%
  dplyr::rename(pi_CHP_Papua_New_Guinea = pi_CHR_PNG) %>%
  dplyr::rename(pi_SAN_Papua_New_Guinea = pi_SAN_PNG) %>%
  dplyr::rename(pi_CHP_Solomon_Islands = pi_CHR_SLM) %>%
  dplyr::rename(pi_SAN_Solomon_Islands = pi_SAN_SLM) %>%
  pivot_longer(cols = starts_with("pi_"), names_to = "pairwise_comb", values_to = "pi_value", names_prefix = "pi_") %>%
  mutate(pop_group = case_when(str_detect(pairwise_comb, "Papua") ~ "Papua New Guinea",
                               str_detect(pairwise_comb, "Solomon") ~ "Solomon Islands", 
                               TRUE ~ "Allopatric"), .before = pi_value)

pi_labels <- pi_df %>% 
  group_by(pairwise_comb) %>%
  dplyr::summarise(
    xpos = max(density(pi_value, na.rm = TRUE)$x[which.max(density(pi_value, na.rm = TRUE)$y)]),
    ypos = max(density(pi_value, na.rm = TRUE)$y))


pi_density_WG <- ggplot(pi_df, aes(x = pi_value, colour = pairwise_comb)) + 
  geom_density(bw = 0.009) +
  coord_cartesian(xlim = c(0,0.3)) +
  labs(x = expression("π"), y = "Density") +
  geom_text_repel(data = pi_labels, aes(x = xpos, y = ypos, label = pairwise_comb), nudge_y = 1) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8))


## plotting everything in a single panel ####

p1 <- (plot_fst / plot_dxy / plot_pi_san / plot_pi_chr / plot_fdM_san / plot_fdM_chr)
p2 <- fst_density_WG / dxy_density_WG / pi_density_WG 
(p1 | p2) + plot_layout(widths = c(2,1))

ggsave("../plots/fig5_1.eps", width = 25, height = 30, units = "cm")
ggsave("../plots/fig5_1.png", width = 25, height = 30, units = "cm")


## alternative with just one fdm
p1 <- (plot_fst / plot_dxy / plot_pi_san / plot_pi_chr / plot_fdM_san)
p2 <- fst_density_WG / dxy_density_WG / pi_density_WG 
(p1 | p2) + plot_layout(widths = c(2,1))

ggsave("../plots/fig5_2.eps", width = 30, height = 30, units = "cm")
ggsave("../plots/fig5_2.pdf", width = 30, height = 30, units = "cm")
ggsave("../plots/fig5_2.svg", width = 30, height = 30, units = "cm")

