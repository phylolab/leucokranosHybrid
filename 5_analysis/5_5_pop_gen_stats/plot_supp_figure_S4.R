###############################|
## 
## Plotting supp. figure S4
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

col.density <- c("#79c99e", "#ff5367", "#0e3b43", "#459aa3", "#FAD4C0", "#b33951", "#ffba49", "#FFA17A", "#7e7f80")
theme_set(theme_light(base_size = 12))

## load data ####
fig_S4_df <- read_csv("../results/all_species_popgen.w50m100.csv.gz") 


## plotting density among populations ####

### fst ####

fst_df <-  fig_S4_df %>%
  select(c(scaffold, start, end, Fst_CHR_FIJI_SAN_PNG, Fst_CHR_FIJI_SAN_SLM, Fst_CHR_FIJI_SAN_AUS, 
           Fst_CHR_SLM_SAN_PNG, Fst_CHR_SLM_SAN_SLM, Fst_CHR_SLM_SAN_AUS, 
           Fst_CHR_PNG_SAN_PNG, Fst_CHR_PNG_SAN_SLM, Fst_CHR_PNG_SAN_AUS)) %>%
  pivot_longer(cols = starts_with("Fst_"), names_to = "pairwise_comb", values_to = "fst_value", names_prefix = "Fst_")

fst_labels <- fst_df %>% 
  group_by(pairwise_comb) %>%
  dplyr::summarise(
    xpos = max(density(fst_value, na.rm = TRUE)$x[which.max(density(fst_value, na.rm = TRUE)$y)]),
    ypos = max(density(fst_value, na.rm = TRUE)$y))

fst_density_WG <- fst_df %>%
  ggplot(aes(fst_value, color = pairwise_comb)) + 
  geom_density(bw = 0.02) +
  coord_cartesian(xlim = c(0.3, 1)) +
  scale_color_manual(values = col.density) +
  labs(x = expression(F[ST]), y = "Density") +
  geom_text_repel(data = fst_labels, aes(x = xpos, y = ypos, label = pairwise_comb), nudge_y = 1) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8))


### dxy ####

dxy_df <- fig_S4_df %>%
  select(c(scaffold, start, end, dxy_CHR_FIJI_SAN_PNG, dxy_CHR_FIJI_SAN_SLM, dxy_CHR_FIJI_SAN_AUS, 
           dxy_CHR_SLM_SAN_PNG, dxy_CHR_SLM_SAN_SLM, dxy_CHR_SLM_SAN_AUS, 
           dxy_CHR_PNG_SAN_PNG, dxy_CHR_PNG_SAN_SLM, dxy_CHR_PNG_SAN_AUS)) %>%
  pivot_longer(cols = starts_with("dxy_"), names_to = "pairwise_comb", values_to = "dxy_value", names_prefix = "dxy_")

dxy_labels <- dxy_df %>% 
  group_by(pairwise_comb) %>%
  dplyr::summarise(
    xpos = max(density(dxy_value, na.rm = TRUE)$x[which.max(density(dxy_value, na.rm = TRUE)$y)]),
    ypos = max(density(dxy_value, na.rm = TRUE)$y))

dxy_density_WG <- dxy_df %>%
  ggplot(aes(dxy_value, color = pairwise_comb)) + 
  geom_density(bw = 0.02) +
  coord_cartesian(xlim = c(0.2, 0.7)) +
  scale_color_manual(values = col.density) +
  labs(x = expression(d[XY]), y = "Density") +
  geom_text_repel(data = dxy_labels, aes(x = xpos, y = ypos, label = pairwise_comb), nudge_y = 1) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8))


## plotting everything in a single panel ####

fst_density_WG / dxy_density_WG
ggsave("../plots/figS4.eps", width = 30, height = 15, units = "cm")
