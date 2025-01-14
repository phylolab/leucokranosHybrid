###################################
###################################
####
#### ADMITXTURE plotting (Alexander et al. 2009)
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
inputprefix <- "results/all_chr_no_lowdp_indv_chryso_sand_leuco.gatk_filtered.no_mono.biall.mac.maxmiss.mindp.maxdp.qual"
inputprefix <- "results/all_chr_no_lowdp_indv_chryso_sand_leuco.gatk_filtered.no_mono.biall.mac.maxmiss.mindp.maxdp.qual.nofiji_aus"

# colour scheme
color.list  <- c("#B33771", "#FC427B", "#BDC581", "#25CCF7", "#F97F51", "#6D214F", "#3B3B98", "#EAB543")
color.fig2  <- c("#cfb27c", "#65433a", "#6c716e", "#8c9560")
color.fig1  <- c("#8b2548", "#44ae9c", "#d5b782")


# read input files -------------------------------------------------------------
qvalue_k2 <- read_table(paste0(inputprefix,".", 2, ".Q"), col_names = F)
qvalue_k3 <- read_table(paste0(inputprefix,".", 3, ".Q"), col_names = F)
qvalue_k4 <- read_table(paste0(inputprefix,".", 4, ".Q"), col_names = F)
qvalue_k5 <- read_table(paste0(inputprefix,".", 5, ".Q"), col_names = F)
qvalue_k6 <- read_table(paste0(inputprefix,".", 6, ".Q"), col_names = F)
qvalue_k7 <- read_table(paste0(inputprefix,".", 7, ".Q"), col_names = F)
crossv    <- read_table(paste0(inputprefix, ".cv.error"), col_names = c("k", "cv_errors")) %>% arrange(cv_errors)
pop       <- read_table(paste0(inputprefix, ".pop"), col_names = c("sample_id", "pop_id"))


# function to plot admixture results -------------------------------------------

plot_admix <- function(pop_file, qvalue_file, K, color_palette, reorder = FALSE){
  #create data frame with population info and qvalues
  qvalue_df <- add_column(pop_file, qvalue_file) %>%
    pivot_longer(cols = starts_with("X"),   #we will convert all columns starting by the name "V" to a single columns
                 names_to = "pop_cluster",  #the name of the columns will be use as variable in a new column names pop_clusters
                 values_to = "qvalue",      #the value within each column will be recoded in a new column called qvalue
                 names_prefix = "X") %>%
    # change order of populations
    mutate(pop_id = factor(pop_id, levels = c("chryso_fiji", "chryso_palau", "chryso_solomon", "chryso_kimbe", "leucokranos", "sanda_kimbe", "sanda_solomon", "sanda_aus")))
  
  # reorder sample_id based on qvalue for pop_cluster = 2 if reorder is TRUE
  if (reorder) {
    # Determine the order based on pop_cluster = 2 qvalue
    sample_order <- qvalue_df %>%
      dplyr::filter(pop_cluster == "2") %>%
      arrange(qvalue) %>%
      pull(sample_id)
    
    # reorder sample_id in qvalue_df based on this custom order
    qvalue_df <- qvalue_df %>%
      mutate(sample_id = factor(sample_id, levels = sample_order))
  }
  # plotting
  ggplot(qvalue_df, aes(fct_relevel(sample_id), qvalue, fill = factor(pop_cluster))) +
    geom_col(color = "gray", size = 0.1, lty = "solid") +
    facet_grid(~pop_id, switch = "x", scales = "free", space = "free") +
    theme_minimal() + 
    labs(x = "Individuals", title = paste("K=", K, sep = ""), y = "Admixture proportion") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete() +
    theme(
      text = element_text(color = "#8d8e8f", family = "Helvetica", size = 15),
      #strip.text = element_blank(), #if we want to remove the name of the population, uncomment this line
      panel.spacing.x = unit(0.5, "lines"),
      panel.border = element_rect(color = "grey25", fill = NA, size = 0.5),
      #axis.text.x = element_blank(), # if we want to remove x axis legend, uncomment this line and comment the next one
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 11),
      axis.ticks.y = element_line(size = 0.2),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.position = "none"
      ) +
    scale_fill_manual(values = color_palette) 
}

# plotting admixture with different k values -----------------------------------

#admix_k2 <- plot_admix(pop, qvalue_k2, 2, color.fig1, reorder = TRUE)
admix_k2 <- plot_admix(pop, qvalue_k2, 2, color.list)
admix_k3 <- plot_admix(pop, qvalue_k3, 3, color.list)
admix_k4 <- plot_admix(pop, qvalue_k4, 4, color.list)
admix_k5 <- plot_admix(pop, qvalue_k5, 5, color.list)
admix_k6 <- plot_admix(pop, qvalue_k6, 6, color.list)
admix_k7 <- plot_admix(pop, qvalue_k7, 7, color.list)

admix_k2 / admix_k3 / admix_k4 / admix_k5 / admix_k6 / admix_k7


# saving plots -----------------------------------------------------------------

ggsave("plots/admix_k2-k7_no_allo.pdf", width = 410, height = 600, units = "mm")
ggsave("plots/admix_k2_no_allo.pdf", width = 280, height = 90, units = "mm")

