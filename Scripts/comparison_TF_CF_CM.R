## TransBiota Project - Transfeminine microbiota analysis
## Obtain abundance figures and core microbiome
## Elaborated by: Jorge Rojas-Vargas


## Load libraries

# Data Wrangling
library(readr)

#

## Read data

setwd("/Users/jorgerv/Documents/transbiota/241019_with_silva138_speciateit2/")
setwd("/Users/jorgerv/Documents/transbiota/241019_with_silva138_speciateit2/R_scripts/")

load("abund_TF_genus.RData")
load("abund_CM_genus.RData")
load("abund_CF_genus.RData")


## --------------------------------------------------------------------------------
##
## Median and prevalence comparing TF, PC, and VC samples
##
## --------------------------------------------------------------------------------

# To obtain data for Tables 2 and S1

# Remove rows named "Unassigned" in all datasets
abund_TF_genus_filtered <- abund_TF_genus[rownames(abund_TF_genus) != "Unassigned", , drop = FALSE]
abund_VC_genus_filtered <- abund_VC_genus[rownames(abund_VC_genus) != "Unassigned", , drop = FALSE]
abund_PC_genus_filtered <- abund_PC_genus[rownames(abund_PC_genus) != "Unassigned", , drop = FALSE]

# Calculate the medians for each genus
medians_TF <- apply(abund_TF_genus_filtered, 1, median)
medians_VC <- apply(abund_VC_genus_filtered, 1, median)
medians_PC <- apply(abund_PC_genus_filtered, 1, median)

# Calculate the percentage of samples with abundance values greater than zero
prevalence_TF <- apply(abund_TF_genus_filtered, 1, function(x) mean(x > 0) * 100)
prevalence_VC <- apply(abund_VC_genus_filtered, 1, function(x) mean(x > 0) * 100)
prevalence_PC <- apply(abund_PC_genus_filtered, 1, function(x) mean(x > 0) * 100)

# Add median and prevalence columns
#abund_TF_genus_filtered$Median <- medians_TF
#abund_TF_genus_filtered$Prevalence <- prevalence_TF

#abund_VC_genus_filtered$Median <- medians_VC
#abund_VC_genus_filtered$Prevalence <- prevalence_VC

#abund_PC_genus_filtered$Median <- medians_PC
#abund_PC_genus_filtered$Prevalence <- prevalence_PC

#abund_TF_taxa$Median <- medians_TF_taxa
#abund_TF_taxa$Prevalence <- prevalence_TF_taxa

# Get the row names with the highest medians
top30_TF <- names(sort(medians_TF, decreasing = TRUE)[1:30])

# Subsample abund_TF, abund_PC, and abund_VC with genera in top30_TF
subsample_TF <- abund_TF_genus_filtered[rownames(abund_TF_genus_filtered) %in% top30_TF, , drop = FALSE]
subsample_PC <- abund_PC_genus_filtered[rownames(abund_PC_genus_filtered) %in% top30_TF, , drop = FALSE]
subsample_VC <- abund_VC_genus_filtered[rownames(abund_VC_genus_filtered) %in% top30_TF, , drop = FALSE]

# Create a function to perform the Wilcoxon test
wilcox_compare <- function(taxa, group1, group2) {
  abund_group1 <- group1[taxa, , drop = FALSE]
  abund_group2 <- group2[taxa, , drop = FALSE]
  p_values <- sapply(taxa, function(taxon) {
    wilcox.test(as.numeric(abund_group1[taxon, ]), as.numeric(abund_group2[taxon, ]))$p.value
  })
  return(p_values)
}

# Compare genera abundances between abund_TF and abund_PC
p_values_TF_PC <- wilcox_compare(top30_TF, abund_TF_genus_filtered, abund_PC_genus_filtered)
names(p_values_TF_PC) <- top30_TF

# Compare taxa abundances between abund_TF and abund_VC
p_values_TF_VC <- wilcox_compare(top30_TF, abund_TF_genus_filtered, abund_VC_genus_filtered)
names(p_values_TF_VC) <- top30_TF

# Create data frames with results
results_TF_PC <- data.frame(Genus = top30_TF, p_value = p_values_TF_PC)
results_TF_PC
results_TF_VC <- data.frame(Genus = top30_TF, p_value = p_values_TF_VC)
results_TF_VC

## Comparing the number of genera among groups

# Count the number of genera (> 0) in each sample for each dataset
genus_counts_TF <- colSums(abund_TF_genus_filtered > 0)
median(genus_counts_TF)
genus_counts_PC <- colSums(abund_PC_genus_filtered > 0)
median(genus_counts_PC)
genus_counts_VC <- colSums(abund_VC_genus_filtered > 0)
median(genus_counts_VC)

# Normality tests (Shapiro-Wilk)
shapiro_TF <- shapiro.test(genus_counts_TF)
shapiro_PC <- shapiro.test(genus_counts_PC)
shapiro_VC <- shapiro.test(genus_counts_VC)

# Display normality test results
cat("Shapiro-Wilk Test:\n")
cat("TF: p-value =", shapiro_TF$p.value, "\n")    # p-value > 0.05 indicates no normal distribution
cat("PC: p-value =", shapiro_PC$p.value, "\n")    # p-value < 0.05 indicates normal distribution
cat("VC: p-value =", shapiro_VC$p.value, "\n\n")  # p-value > 0.05 indicates no normal distribution

# Compare TF vs PC using the Wilcoxon test
test_TF_PC <- wilcox.test(genus_counts_TF, genus_counts_PC, alternative = "two.sided")
test_TF_PC 

# Compare TF vs VC using the Wilcoxon test
test_TF_VC <- wilcox.test(genus_counts_TF, genus_counts_VC, alternative = "two.sided")
test_TF_VC
