## TransBiota Project - Transfeminine microbiota analysis
## Obtain abundance figures and the core microbiome
## Elaborated by: Jorge Rojas-Vargas


## Load libraries

# Explore microbiome
library(phyloseq)
library(microbiome)

# Plot tools
library(RColorBrewer)
library(viridis)
library(svglite)

## Read data

load("abund_TF_genus.RData")
load("metadata.RData")
metadata <- as.data.frame(metadata)


## --------------------------------------------------------------------------------
##
## Core microbiome of all TF samples
##
## --------------------------------------------------------------------------------


# Eliminate "Unassigned" taxa from genus abundance table in TF samples
abund_TF_genus_filtered <- abund_TF_genus[row.names(abund_TF_genus) != "Unassigned", ]

# ps object of abundance table
taxa_table <- as.matrix(row.names(abund_TF_genus_filtered))
colnames(taxa_table) <- "Genus"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(abund_TF_genus_filtered, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
sampledata <- sample_data(metadata)
rownames(sampledata) <- metadata$UID

ps_abund_genus_filt = phyloseq(OTU_taxa, TAX_taxa, sampledata)
ps_abund_genus_filt <- prune_taxa(taxa_sums(ps_abund_genus_filt)>0, ps_abund_genus_filt)
ps_abund_genus_filt <- microbiome::transform(ps_abund_genus_filt, "compositional")

core_TF_genus <- core_members(ps_abund_genus_filt, detection = 1/100, prevalence = 50/100, include.lowest = TRUE)
print(core_TF_genus)
ps.core_TF_genus <- core(ps_abund_genus_filt, detection = 1/100, prevalence = 50/100, include.lowest = TRUE)

prevalences <- seq(5/100, 100/100, 5/100)
detections <- round(10^seq(log10(0.001), log10(.5), length = 10), 3)

p1 <- plot_core(ps.core_TF_genus,
                plot.type = "heatmap",
                colours = rev(brewer.pal(5, "RdGy")),
                prevalences = prevalences,
                detections = detections, 
                min.prevalence = 0.01) +
  xlab("Detection Threshold (Relative Abundance)") +
  theme(axis.text.y = element_text(size = 10, face = "italic"),
        axis.text.x.bottom = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

p1 <- p1 + theme_bw() + ylab("Genus")
p1


print(p1 + scale_fill_viridis())

# svglite("Core_microbiome_of_Genus.svg", width = 6.5, height = 3)
(p1 + scale_fill_viridis())
# dev.off()
