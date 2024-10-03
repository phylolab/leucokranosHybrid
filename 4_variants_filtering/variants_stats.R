##### The script is used to examine the statistics of variants #####

# Load the required package
# install.packages("tidyverse") # Uncomment this line if the package hasn't been installed
# install.packages("ggplot2") # Uncomment this line if the package hasn't been installed
library(tidyverse)
library(ggplot2)

# Set the variables
#!! Change for each dataset
inputprefix <- "AllSamples_SNP.SG.Filtered"

# Create a folder to output the stats plots
dir.create(paste0(inputprefix,".Stats"), recursive = T)

############################
# Variant based statistics #
############################
# 1. Variant quality (Phred encoded site quality)
# Read the file with site quality
variant_qual <- read.delim(paste0(inputprefix,".lqual"),
                          sep = "\t",
                          col.names = c("Chr","POS","QUAL"),
                          header = T)

fig <- ggplot(variant_qual, aes(QUAL)) +
        geom_density(fill = "lightblue", colour = "black", alpha = 0.3) +
        ggtitle("Distribution of Variant Quality") +
        xlab("Variant Quality (Phred score)") +
        ylab("Density") +
        labs(caption = paste("Dataset:",inputprefix))
fig + theme_light() + xlim(0,500) 
    
ggsave(paste0(inputprefix,".Stats/DistributionVariantQUAL.png"), plot = last_plot())

# 2. Variant mean depth
# Read the file with variant mean depth
variant_depth <- read.delim(paste0(inputprefix,".ldepth.mean"),
                            sep = "\t",
                            col.names = c("Chr","POS","Mean_Depth","Var_Depth"),
                            header = T)

# First get an idea about how data distributes 
summary(variant_depth$Mean_Depth)

# Plot the distribution with limited x range inferred from the summary
fig <- ggplot(variant_depth, aes(Mean_Depth)) +
      geom_density(fill = "lightblue", colour = "black", alpha = 0.3) +
      geom_vline(xintercept = 6, colour = "red", linetype = "dashed") +
      ggtitle("Distribution of Variant Mean Depth") +
      xlab("Variant Mean Depth") +
      ylab("Density") +
      labs(caption = paste("Dataset:",inputprefix))
fig + theme_light() + xlim(0,25) 

ggsave(paste0(inputprefix,".Stats/DistributionVariantMeanDepth.png"), plot = last_plot())


# 3. Variant missingness
# Read the file with variant missingness
variant_missing <- read.delim(paste0(inputprefix,".lmiss"),
                              sep = "\t",
                              col.names = c("Chr","POS","nChr","nFiltered","nMiss","fMiss"),
                              header = T)

# First get an idea about how data distributes
summary(variant_missing$fMiss)
  #! fMiss stands for the fraction of samples with missing data for each site

# Plot the distribution
fig <- ggplot(variant_missing, aes(fMiss))  +
  geom_density(fill = "lightblue", colour = "black", alpha = 0.3) +
  ggtitle("Distribution of Variant Missingness") +
  xlab("Fraction of Samples with Missing Genotypes") +
  ylab("Density") +
  labs(caption = paste("Dataset:",inputprefix))
fig + theme_light() + xlim(0,0.05)

ggsave(paste0(inputprefix,".Stats/DistributionVariantMissing.png"), plot = last_plot())


# 4. Minor allele frequency
# Read the file with minor-allele frequency (MAF)
variant_MAF <- read.delim(paste0(inputprefix, ".frq"),
                          sep = "\t",
                          col.names = c("Chr","POS","nAlleles","nChr","A1","A2"),
                          header = T)

# Find minor allele frequency
variant_MAF$MAF <- variant_MAF %>% select(A1,A2) %>% apply(1, function(afrq) min(afrq))

# First get an idea about how data distributes
summary(variant_MAF$MAF)

# Plot the distribution of MAF
fig <- ggplot(variant_MAF, aes(MAF))  +
  geom_density(fill = "lightblue", colour = "black", alpha = 0.3) +
  ggtitle("Distribution of Minor-Allele Frequency") +
  xlab("Minor-Allele Frequency") +
  ylab("Density") +
  labs(caption = paste("Dataset:",inputprefix))
fig + theme_light() + xlim(0,0.5)

ggsave(paste0(inputprefix,".Stats/DistributionMAF.png"), plot = last_plot())



###############################
# Individual based statistics #
###############################
# 1. Sample mean depth
# Read the file with mean depth per sample
sample_depth <- read.delim(paste0(inputprefix, ".idepth"),
                           sep = "\t",
                           col.names = c("Sample","nSites","Depth"),
                           header = T)


# Plot the distribution of mean depth
fig <- ggplot(sample_depth, aes(Depth)) +
    geom_histogram(fill = "navy", colour = "black", alpha = 0.3, bins = 50) +
    geom_vline(xintercept = 6, colour = "red", linetype = "dashed") +
    ggtitle("Distribution of Individual Mean Depth") +
    xlab("Individual Mean Depth") +
    ylab("Count") +
    labs(caption = paste("Dataset:",inputprefix))
fig + theme_light() 

ggsave(paste0(inputprefix,".Stats/DistributionIndMeanDepth.png"),plot = last_plot())

# 2. Sample missingness
# Read the file with sample missingness
sample_missing <- read.delim(paste0(inputprefix, ".imiss"),
                             sep = "\t",
                             col.names = c("Sample","nData","nFiltered","nMiss","fMiss"),
                             header = T)

# First get an idea about how data distributes
summary(sample_missing$fMiss)

# Plot the distribution
fig <- ggplot(sample_missing, aes(fMiss))  +
  geom_histogram(fill = "navy", colour = "black", alpha = 0.3) +
  ggtitle("Distribution of Individual Missingness") +
  xlab("Fraction of Sites with Missing Genotypes") +
  ylab("Count") +
  labs(caption = paste("Dataset:",inputprefix))
fig + theme_light() + xlim(0,0.01)

ggsave(paste0(inputprefix,".Stats/DistributionIndMissing.png"), plot = last_plot())

# 3. Heterozygosity and inbreeding coefficient per individual
# Read the file heterozygosity and inbreeding coefficient
sample_het <- read.delim(paste0(inputprefix, ".het"),
                         sep = "\t",
                         col.names = c("Sample","HO","HE","nSites","F"),
                         header = T)

# Plot the distribution
fig <- ggplot(sample_het, aes(F))  +
  geom_histogram(fill = "navy", colour = "black", alpha = 0.3, bins = 50) +
  ggtitle("Distribution of Individual Inbreeding Coefficient") +
  xlab("Inbreeding Coefficient") +
  ylab("Count") +
  labs(caption = paste("Dataset:",inputprefix))
fig + theme_light() 

ggsave(paste0(inputprefix,".Stats/DistributionIndInbreedCoef.png"), plot = last_plot())




