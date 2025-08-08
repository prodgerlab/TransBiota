## TransBiota Project - Transfeminine microbiota analysis
## Calculate alpha and beta diversity
## Elaborated by: Jorge Rojas-Vargas


## Load libraries

#Explore microbiome
library(phyloseq)
library(vegan)

#Plot tools
library(svglite)

#Statistical tools
library(RVAideMemoire)


## Read data

#It is important to use the untouched count table here as we're interested in rare taxa.
load("count_TF_genus.RData")
load("count_TF_taxa.RData")
load("count_PC_genus.RData")
load("count_PC_taxa.RData")
load("count_VC_genus.RData")
load("count_VC_taxa.RData")

#Load phyloseq files with CLR transformation (via clr_c function of the Tjazi package)
load("ps_count_TF_genus.clr.RData")



## --------------------------------------------------------------------------------
##
## Alpha diversity comparing TF, CM, and CF samples
##
## --------------------------------------------------------------------------------


#Compute alpha diversity using the vegan library to evaluate normality

# Define groups and count tables
groups <- c("TF", "PC", "VC")
count_tables <- list(TF = count_TF_genus, PC = count_PC_genus, VC = count_VC_genus)
#count_tables <- list(TF = count_TF_taxa, PC = count_PC_taxa, VC = count_VC_taxa)

# Compute diversity indices and observed taxa
diversity_indices <- lapply(count_tables, function(ct) {
  list(
    shannon = vegan::diversity(t(ct), index = "shannon"),
    simpson = vegan::diversity(t(ct), index = "simpson"),
    observed = rowSums(t(ct) > 0)
  )
})

# Run Shapiro-Wilk test for each index in each group
shapiro_results <- lapply(diversity_indices, function(indices) {
  sapply(indices, shapiro.test, simplify = FALSE)
})

# Combine all indices into a data frame
diversity_data <- do.call(rbind, lapply(groups, function(group) {
  data.frame(
    group = group,
    shannon = diversity_indices[[group]]$shannon,
    simpson = diversity_indices[[group]]$simpson,
    observed = diversity_indices[[group]]$observed
  )
}))

# Pairwise Wilcoxon Tests
wilcox_results <- lapply(c("shannon", "simpson", "observed"), function(index) {
  pairwise.wilcox.test(
    diversity_data[[index]],
    diversity_data$group,
    p.adjust.method = "fdr"
  )
})

# Print results
names(wilcox_results) <- c("Shannon", "Simpson", "Observed")
lapply(names(wilcox_results), function(name) {
  cat(paste("\n", name, "Index Pairwise Wilcoxon Test:\n", sep = ""))
  print(wilcox_results[[name]])
})


# Calculate medians for Shannon Index
median_shannon_TF <- median(diversity_indices$TF$shannon)
median_shannon_PC <- median(diversity_indices$PC$shannon)
median_shannon_VC <- median(diversity_indices$VC$shannon)

# Print results
cat("Median Shannon Index:\n")
cat("TF:", median_shannon_TF, "\n")
cat("PC:", median_shannon_PC, "\n")
cat("VC:", median_shannon_VC, "\n")






## --------------------------------------------------------------------------------
##
## Beta diversity of TF samples
##
## --------------------------------------------------------------------------------

### We choose the Aitchison distance using CLR transformed data at the GENUS level


#Generate distance matrix using Aitchison distances
dist_matrix <- phyloseq::distance(ps_count_TF_genus.clr, method = "euclidean")


## ------------------------------------------
## Dispersion test with BCs
## ------------------------------------------

#Calculate the beta dispersion
dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(ps_count_TF_genus.clr)$BC)

#Check the results of the dispersion test
dispr

#PERMANOVA test to dispersion
#Set a seed for reproducibility
set.seed(123)
vegan::adonis2(dist(dispr$distances) ~dispr$group) # Global p-value 0.191

#Pairwise PERMANOVA test with fdr correction of dispersion
set.seed(123)
pairwise_test <- vegan::permutest(dispr, pairwise = TRUE)
pairwise_pvalues <- pairwise_test$pairwise$permuted 
pairwise_pvalues
pairwise_pvalues_fdr <- p.adjust(pairwise_pvalues, method = "fdr")
pairwise_pvalues_fdr 
# All q-values > 0.2

# PERMANOVA test for microbiota beta diversity (community structure)
set.seed(123)
vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps_count_TF_genus.clr)$BC)
# global p-value 0.001

#Pairwise PERMANOVA test with fdr correction for microbiota beta diversity (community structure)
set.seed(123)
pairwise.perm.manova(dist_matrix, phyloseq::sample_data(ps_count_TF_genus.clr)$BC, nperm = 10000)
# BC1-BC2 q-value 0.030
# BC1-BC4 q-value 0.026
# BC2-BC4 q-value 0.030


#PCA and dispersion plot

#svglite("disp_pca_BC_genus.svg", width=10, height=7)
plot(dispr, main = "", sub = "",
     xlab = "PCA 1 (16.26%)", ylab = "PCA 2 (12.02%)",
     col = c(BC_1 = "#fb9d10", 
             BC_2 = "#ff0000",
             BC_3 = "#006efe", 
             BC_4= "black"
     ), 
     cex = 1.5)
#dev.off()




## ------------------------------------------
## Dispersion test with Symptoms
## ------------------------------------------

dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(ps_count_TF_genus.clr)$Symptom)
#PCA1 16.55%, PCA2 12.02%

#Check the results of the dispersion test
dispr

# PERMANOVA test for dispersion
#Set a seed for reproducibility
set.seed(123)
vegan::adonis2(dist(dispr$distances) ~dispr$group)
# Global p-value 0.001

#Pairwise PERMANOVA test with fdr correction for dispersion
set.seed(123)
pairwise_test <- vegan::permutest(dispr, pairwise = TRUE)
pairwise_pvalues <- pairwise_test$pairwise$permuted 
pairwise_pvalues_fdr <- p.adjust(pairwise_pvalues, method = "fdr")
pairwise_pvalues_fdr
# Discharge-Itch_Malod q-value 0.021
# Malodour-None q-value 0.021

# Pairwise comparisons of microbiota beta diversity (community structure)
set.seed(123)
vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps_count_TF_genus.clr)$Symptom)
# global p-value 0.006

#Pairwise PERMANOVA test with fdr correction for microbiota beta diversity (community structure)
set.seed(123)
pairwise.perm.manova(dist_matrix, phyloseq::sample_data(ps_count_TF_genus.clr)$Symptom, nperm = 10000)
# Malodour-Bleeding q-value 0.0098
# Malodour-Itch q-value 0.0147
# Malodour-None q-value 0.0180


#PCA and dispersion plot
svglite("disp_pca_TF_genus.svg", width=10, height=7)
plot(dispr, main = "", sub = "",
     xlab = "PCA 1 (16.55%)", ylab = "PCA 2 (12.02%)",
     col = c(Bleeding = "#fb9d10", 
             Bleed_Disch = "#9f69ff", 
             Bleed_Malod = "#FF1FB4",
             Bleed_Pain = "#ff0000",
             Discharge = "#026b6b", 
             Itch = "#006efe", 
             Itch_Malod = "#00e6e6",
             Malod_Disch = "#7ca386",
             Malodour = "#f5d104", 
             None = "black", 
             Pain = "#666666"
     ), 
     cex = 1.5)
dev.off()

#vglite("disp_boxplot_TF.svg", width=10, height=7)
boxplot(dispr, xlab = "", 
        col = c("#fe9929", "#8c6bb1", "#2a9d8f"))
#dev.off()
