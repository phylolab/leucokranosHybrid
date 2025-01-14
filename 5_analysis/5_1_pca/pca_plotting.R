###################################
###################################
####
#### PCA (with plink output)
#### on pruned dataset
####
#### Author: Sarah Schmid
####
###################################
###################################

# load the required packages ---------------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(here)


# set the variables ------------------------------------------------------------
inputprefix <- "../results/all_chr_no_lowdp_indv_chryso_sand_leuco.gatk_filtered.no_mono.biall.mac.maxmiss.mindp.maxdp.qual"
#inputprefix <- "../results/all_chr_no_lowdp_indv_all.gatk_filtered.no_mono.biall.mac.maxmiss.mindp.maxdp.qual"

color.list <- c("#B33771", "#FC427B", "#BDC581", "#25CCF7", "#F97F51", "#6D214F", "#3B3B98", "#EAB543")
color.fig2 <- c("#cfb27c", "#65433a", "#6c716e", "#8c9560")
color.fig1 <- c("#44ae9c", "#d5b782", "#8b2548")


# load data and prepare dataframe for plotting ---------------------------------
pca <- read_table(paste0(inputprefix,".eigenvec"), col_names = F)
eigenval <- scan(paste0(inputprefix, ".eigenval"))

# remove extra ID column
pca <- pca[,-1]

# set column names
colnames(pca)[1] <- "SampleID"
colnames(pca)[2:ncol(pca)] <- paste0("PC",1:(ncol(pca)-1))

# percentage variance explained
pve <- data.frame(PC = 1:(ncol(pca)-1), pve = eigenval/sum(eigenval)*100)

# add species information as we also want to see whether variants cluster by species
sp <- rep(NA, length(pca$SampleID))
sp[grep("CH", pca$SampleID)] <- "A.chrysopterus"
sp[grep("LU", pca$SampleID)] <- "A.leucokranos"
sp[grep("SA", pca$SampleID)] <- "A.sandaracinos"
sp[grep("GB", pca$SampleID)] <- "A.clarkii"
pca <- as_tibble(data.frame(pca, sp))


# plotting data ----------------------------------------------------------------

## plot PCA for hybrids + parents ----------------------------------------------
fig <- pca[which(!pca$sp=="A.clarkii"),] %>% 
  ggplot(aes(x = PC1, y = PC2, colour = sp, label = SampleID)) +
  geom_point(size = 5) + 
  scale_colour_manual(values = c("#8b2548", "#d5b782", "#44ae9c")) +
  xlab(paste0("PC1 (",round(pve[1,2], digits = 2),"%)")) +
  ylab(paste0("PC2 (",round(pve[2,2], digits = 2),"%)")) +
  geom_label_repel(aes(label = SampleID), max.overlaps = 50, alpha = 0.6, size = 2) +
  theme_light() + 
  theme(
    legend.position = "none",
    text = element_text(color = "#8d8e8f", family = "Helvetica", size = 15),
    axis.text = element_text(color = "#8d8e8f", size = 16),
    axis.title.x.bottom = element_text(vjust = 14, hjust = 0.98, size = 16, color = "#8d8e8f"),
    axis.title.y.left = element_text(vjust = -17, hjust = 0.98, size = 16, color = "#8d8e8f"),
    panel.grid.minor = element_blank()
  )

ggsave("../plots/pca_parents_hybrid.pdf", width = 180, height = 180, units = "mm")
ggsave("../plots/pca_parents_hybrid_labels.pdf", width = 180, height = 180, units = "mm")


### plot PCA for hybrids + parents (PCA1 + 3) ----------------------------------------------
fig_pc1_3 <- pca[which(!pca$sp=="A.clarkii"),] %>% 
  ggplot(aes(x = PC1, y = PC3, colour = sp, label = SampleID)) +
  geom_point(size = 5) + 
  scale_colour_manual(values = c("#8b2548", "#d5b782", "#44ae9c")) +
  xlab(paste0("PC1 (",round(pve[1,2], digits = 2),"%)")) +
  ylab(paste0("PC3 (",round(pve[3,2], digits = 2),"%)")) +
  geom_label_repel(aes(label = SampleID), max.overlaps = 50, alpha = 0.6, size = 2) +
  theme_light() + 
  theme(
    legend.position = "none",
    text = element_text(color = "#8d8e8f", family = "Helvetica", size = 15),
    axis.text = element_text(color = "#8d8e8f", size = 16),
    axis.title.x.bottom = element_text(vjust = 13, hjust = 0.02, size = 16, color = "#8d8e8f"),
    axis.title.y.left = element_text(vjust = -16, hjust = 0.50, size = 16, color = "#8d8e8f"),
    panel.grid.minor = element_blank()
  )

ggsave("../plots/pca_parents_hybrid_pc1_3.pdf", width = 180, height = 180, units = "mm")
ggsave("../plots/pca_parents_hybrid_pc1_3_labels.pdf", width = 180, height = 180, units = "mm")


### plot PCA for hybrids + parents (PCA1 + 4) ----------------------------------------------
fig_pc1_4 <- pca[which(!pca$sp=="A.clarkii"),] %>% 
  ggplot(aes(x = PC1, y = PC4, colour = sp, label = SampleID)) +
  geom_point(size = 5) + 
  scale_colour_manual(values = c("#8b2548", "#d5b782", "#44ae9c")) +
  xlab(paste0("PC1 (",round(pve[1,2], digits = 2),"%)")) +
  ylab(paste0("PC4 (",round(pve[4,2], digits = 2),"%)")) +
  geom_label_repel(aes(label = SampleID), max.overlaps = 50, alpha = 0.6, size = 2) +
  theme_light() + 
  theme(
    legend.position = "none",
    text = element_text(color = "#8d8e8f", family = "Helvetica", size = 15),
    axis.text = element_text(color = "#8d8e8f", size = 16),
    axis.title.x.bottom = element_text(vjust = 13, hjust = 0.02, size = 16, color = "#8d8e8f"),
    axis.title.y.left = element_text(vjust = -15, hjust = 0.50, size = 16, color = "#8d8e8f"),
    panel.grid.minor = element_blank()
  )

ggsave("../plots/pca_parents_hybrid_pc1_4.pdf", width = 180, height = 180, units = "mm")
ggsave("../plots/pca_parents_hybrid_pc1_4_labels.pdf", width = 180, height = 180, units = "mm")

### combine both plots together for supp. mat.

fig_pc1_3 + fig_pc1_4 + plot_annotation(tag_levels = 'A')
ggsave("../plots/fig_s1_pca1_3_pca1_4_with_labels.pdf", width = 360, height = 180, units = "mm")


## plot PCA for leucokranos ----------------------------------------------------

# select the right dataset
inputprefix <- "../results/all_chr_no_lowdp_indv_leuco_filtered.gatk_filtered.no_mono.biall.mac.maxmiss.mindp.maxdp.qual"

# load data and prepare dataframe for plotting ---------------------------------
pca <- read_table(paste0(inputprefix,".eigenvec"), col_names = F)
eigenval <- scan(paste0(inputprefix, ".eigenval"))

# remove extra ID column
pca <- pca[,-1]

# set column names
colnames(pca)[1] <- "SampleID"
colnames(pca)[2:ncol(pca)] <- paste0("PC",1:(ncol(pca)-1))

# percentage variance explained
pve <- data.frame(PC = 1:(ncol(pca)-1), pve = eigenval/sum(eigenval)*100)

# add species information as we also want to see whether variants cluster by species
sp <- rep(NA, length(pca$SampleID))
sp[grep("LU", pca$SampleID)] <- "A.leucokranos"
pca <- as_tibble(data.frame(pca, sp))


fig_pca_leuco <- pca[which(pca$sp=="A.leucokranos"),] %>% 
  ggplot(aes(x = PC1, y = PC2, colour = sp, label = SampleID)) +
  geom_point(size = 5) + 
  scale_colour_manual(values = "#d5b782") +
  xlab(paste0("PC1 (",round(pve[1,2], digits = 2),"%)")) +
  ylab(paste0("PC2 (",round(pve[2,2], digits = 2),"%)")) +
  #geom_label_repel(aes(label = SampleID), max.overlaps = 50, alpha = 0.6, size = 2) +
  geom_text_repel(aes(label = SampleID),
                   size = 5,
                   box.padding   = 0.35,
                   point.padding = 0.8,
                   segment.color = "grey50",
                   max.overlaps = 35,
                   show.legend = F) +
  theme_light() + 
  theme(
    legend.position = "none",
    text = element_text(color = "#8d8e8f", family = "Helvetica", size = 15),
    axis.text = element_text(color = "#8d8e8f", size = 16),
    axis.title.x.bottom = element_text(vjust = 13, hjust = 0.01, size = 16, color = "#8d8e8f"),
    axis.title.y.left = element_text(vjust = -16, hjust = 0.99, size = 16, color = "#8d8e8f"),
    panel.grid.minor = element_blank()
  )

ggsave("../plots/pca_hybrids_pc1_2.eps", width = 180, height = 180, units = "mm")



