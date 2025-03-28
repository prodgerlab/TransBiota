## TransBiota Project - Transfeminine microbiota analysis
## To generate the stack bar plots of the paper
## Elaborated by: Jorge Rojas-Vargas


## Load libraries

# Explore microbiome
library(ggplot2)
library(svglite)
library(dplyr)
# Graphical tools
library(reshape2)
library(patchwork)

## Read data
# TF abundance
load("abund_TF_genus.RData")
# rCF abundance
load("abund_CF_genus.RData")
# CM abundance
load("abund_CM_genus.RData")
# Metadata
load("metadata.RData")
# Samples ID
load("samples_id.RData")


## --------------------------------------------------------------------------------
##
## Stackbar plots Fig 2 and Fig S2
##
## --------------------------------------------------------------------------------

# Change the name of S5-A14a
rownames(abund_TF_genus)[rownames(abund_TF_genus) == "S5-A14a"] <- "S5_A14a"
rownames(abund_VC_genus)[rownames(abund_VC_genus) == "S5-A14a"] <- "S5_A14a"
rownames(abund_PC_genus)[rownames(abund_PC_genus) == "S5-A14a"] <- "S5_A14a"


## For baseline comparisons
# Filter rows that have "W1" in the week column
metadata_filtered <- metadata[metadata$Week == "W1", ]
# # Add TMI032_W2 and TMI038_W2 to metadata_filtered
# Get unique identifiers from metadata
valid_samples_baseline <- metadata_filtered$UID 
# Add TMI032_W2 and TMI038_W2 to valid_samples
valid_samples_baseline <- c(valid_samples_baseline, "TMI032_W2", "TMI038_W2") # Generate a full count table including only the baseline (week 1) for the TF samples
valid_samples_PC <- samples_id %>% dplyr::filter(Type == "PC") %>% pull(Name)  
valid_samples_VC <- samples_id %>% dplyr::filter(Type == "VC") %>% pull(Name) 
# Combine the three lists into one and get only unique elements
all_valid_samples_baseline <- unique(c(valid_samples_PC, valid_samples_VC, valid_samples_baseline))


### Top taxa according to the median values

# Selection of method to determine the top taxa
method <- "median"    # using median values as selection criteria

## Selection of taxa level
# Eliminate "Unassigned" at Genus level
abund_TF_filt <- abund_TF_genus[rownames(abund_TF_genus) != "Unassigned", ]
abund_PC_filt <- abund_PC_genus[rownames(abund_PC_genus) != "Unassigned", ]
abund_VC_filt <- abund_VC_genus[rownames(abund_VC_genus) != "Unassigned", ]

# Eliminate genus "Bacteria" at Taxa level
# abund_TF_filt <- abund_TF_genus[rownames(abund_TF_genus) != "Bacteria", ]
# abund_PC_filt <- abund_PC_genus[rownames(abund_PC_genus) != "Bacteria", ]
# abund_VC_filt <- abund_VC_genus[rownames(abund_VC_genus) != "Bacteria", ]

# Calculate the value according
TF_filt <- apply(abund_TF_filt, 1, method)
PC_filt <- apply(abund_PC_filt, 1, method)
VC_filt <- apply(abund_VC_filt, 1, method)


# Get the row names with the highest medians
top_TF <- names(sort(TF_filt, decreasing = TRUE)[1:30])
top_PC <- names(sort(PC_filt, decreasing = TRUE)[1:30])
top_VC <- names(sort(VC_filt, decreasing = TRUE)[1:30])

# Combine the taxonomies and get the unique ones
top_taxa <- unique(c(top_TF, top_PC, top_VC))

# Subsamples of each abundance table using TF baseline
abund_TF_baseline <- abund_TF_filt[, which(colnames(abund_TF_filt) %in% valid_samples_baseline)] # TF baseline

top_abund_TF <- abund_TF_baseline[rownames(abund_TF_baseline) %in% top_taxa, ] # if want to plot just TF baseline
top_abund_TF <- abund_TF_filt[rownames(abund_TF_filt) %in% top_taxa, ] # if want to plot all TF samples

top_abund_PC <- abund_PC_filt[rownames(abund_PC_filt) %in% top_taxa, ]
top_abund_VC <- abund_VC_filt[rownames(abund_VC_filt) %in% top_taxa, ]

# Add a row called "Other" in each subsample
other_TF <- 1 - colSums(top_abund_TF)
other_PC <- 1 - colSums(top_abund_PC)
other_VC <- 1 - colSums(top_abund_VC)
Chryseobacterium_TF <- colSums(top_abund_TF) - colSums(top_abund_TF)

top_abund_TF <- rbind(top_abund_TF, Chryseobacterium = Chryseobacterium_TF, Other = other_TF)
top_abund_PC <- rbind(top_abund_PC, Other = other_PC)
top_abund_VC <- rbind(top_abund_VC, Other = other_VC)

# Melt the subsamples to convert from wide to long format
dat_TF <- melt(as.matrix(top_abund_TF))
dat_PC <- melt(as.matrix(top_abund_PC))
dat_VC <- melt(as.matrix(top_abund_VC))

# Rename the columns
colnames(dat_TF) <- c("Genus", "Sample", "Abundance")
colnames(dat_PC) <- c("Genus", "Sample", "Abundance")
colnames(dat_VC) <- c("Genus", "Sample", "Abundance")


# Extract the medians corresponding to the taxonomies in top_taxa
median_top_taxa <- TF_filt[top_taxa]
# Order top_taxa by the medians, from highest to lowest
top_taxa_order <- names(sort(median_top_taxa, decreasing = TRUE))
# Add missing taxa and "Other" to the end of top_taxa_order
top_taxa_order <- c(top_taxa_order, "Chryseobacterium", "Other")
# Move "Prevotella" to the first place and "Peptoniphilus" to the second
top_taxa_order <- c("Prevotella", "Peptoniphilus", setdiff(top_taxa_order, c("Prevotella", "Peptoniphilus")))
# Show the final result of top_taxa_order
top_taxa_order

# Convert the Genus column to a factor with levels in top_taxa_order
dat_TF$Genus <- factor(dat_TF$Genus, levels = top_taxa_order)
dat_PC$Genus <- factor(dat_PC$Genus, levels = top_taxa_order)
dat_VC$Genus <- factor(dat_VC$Genus, levels = top_taxa_order)

# Order the dataframes according to the factor levels of Genus
dat_TF <- dat_TF[order(dat_TF$Genus), ]
dat_PC <- dat_PC[order(dat_PC$Genus), ]
dat_VC <- dat_VC[order(dat_VC$Genus), ]

# Organize the data to be in ascending order by abundance of X taxa
X = "Prevotella"
dat_TF_X <- subset(dat_TF, Genus == X)
dat_TF_X <- dat_TF_X[order(dat_TF_X$Abundance), ]
ordered_Sample <- dat_TF_X$Sample
dat_TF$Sample <- factor(dat_TF$Sample, levels = ordered_Sample)
#dat_tF$Genus <- factor(dat_tF$Genus, levels = c(X, setdiff(levels(dat_tF$Genus), X)))

dat_PC_X <- subset(dat_PC, Genus == X)
dat_PC_X <- dat_PC_X[order(dat_PC_X$Abundance), ]
ordered_Sample <- dat_PC_X$Sample
dat_PC$Sample <- factor(dat_PC$Sample, levels = ordered_Sample)

X = "Lactobacillus"
dat_VC_X <- subset(dat_VC, Genus == X)
dat_VC_X <- dat_VC_X[order(dat_VC_X$Abundance), ]
ordered_Sample <- dat_VC_X$Sample
dat_VC$Sample <- factor(dat_VC$Sample, levels = ordered_Sample)
dat_VC$Genus <- factor(dat_VC$Genus, levels = c("Lactobacillus", setdiff(levels(dat_VC$Genus), "Lactobacillus")))


# Define color palette for taxa
my_colors <- c(
  Acinetobacter = "#b78bf0",
  Actinotignum = "#355c7d",
  Anaerococcus = "#87cefa",
  Arcanobacterium = "#6209d6",
  Atopobium = "#9f69ff",
  Bacteroides = "#782144",
  Berryella = "#ff1292",
  Campylobacter = "#34eb9b",
  Chryseobacterium = "#f795ff",
  Corynebacterium = "#ffff00",
  Cutibacterium = "#a02c5a",
  Dialister = "#997dbd",
  Ezakiella = "#026b6b",
  Fannyhessea = "#ffffb5",
  Fastidiosipila = "#c5be63",
  Fenollaria = "#f5d104",
  Finegoldia = "#800080",
  Fusobacterium = "#CD853F",
  Gardnerella = "#784421",
  Howardella = "#1FFF87",
  Hoylesella = "#0989d7",
  HT002 = "#fb9d10",
  Lactobacillus = "#ff0000",
  Lawsonella = "#1FFF9E",
  Megasphaera = "#0000ff",
  Mobiluncus = "#f08080",
  Murdochiella = "#ca03fc",
  Negativicoccus = "#00e6e6",
  Parvimonas = "#00ffff",
  Peptococcus = "#5e4fa2",
  Peptoniphilus = "#cccc00",
  Peptostreptococcus = "#DEDB10",
  Porphyromonas = "#f16e43",
  Prevotella = "#6bbb88",
  Propionimicrobium = "#8F00FF",
  Roseateles = "#006efe",
  S5_A14a = "#ac939d",  # Note: changed "-" to "_" for valid name
  Schaalia = "#a2ffb3",
  Sneathia = "#d1e8eb",
  Staphylococcus = "#7ca386",
  Streptococcus = "#ffc0cb",
  Ureaplasma = "#290B00",
  Varibaculum = "#094717",
  Veillonella = "#499996",
  W5053 = "#FF1FB4",
  Winkia = "#ffdfa2",
  Other = "#666666"
)



#Create a correspondence between Genus and color
# genus_colors <- setNames(my_colors, top_taxa_order)

#Add a column indicating the origin of each dataframe
dat_TF$Dataset <- "dat_TF"
dat_PC$Dataset <- "dat_PC"
dat_VC$Dataset <- "dat_VC"

# Combine the dataframes into one
combined_data <- bind_rows(dat_TF, dat_PC, dat_VC)

# Create individual plots with manually assigned colors
plot_TF <- ggplot(dat_TF, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") + 
  xlab("") + 
  ylab("Relative Abundance") + 
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.position = "none",
    #legend.text = element_text(size = 14), legend.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  ) + 
  scale_fill_manual(values = my_colors) +
  facet_wrap(~ Dataset)

plot_PC <- ggplot(dat_PC, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") + 
  xlab("") + 
  ylab("") + 
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) + 
  scale_fill_manual(values = my_colors) +
  facet_wrap(~ Dataset)

plot_VC <- ggplot(dat_VC, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") + 
  xlab("Sample") + 
  ylab("") + 
  theme_bw() +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    legend.position = "none"
  ) + 
  scale_fill_manual(values = my_colors) +
  facet_wrap(~ Dataset)

(plot_TF | plot_PC | plot_VC)

#svglite("Fig_2.svg", width=20, height=10) # Showing only TF baseline
(plot_TF | plot_PC | plot_VC)
#dev.off()

#svglite("Supp_Fig_2.svg", width=20, height=10) # Showing all TF samples
plot_TF
#dev.off()


#svglite("legend_fig_2.svg", width=10, height=8)
ggplot(dat_TF, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") + xlab("Sample") + ylab("Abundance") + theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12), axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"), 
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 16, face = "bold")
  ) +  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) + 
  scale_fill_manual(
    values = my_colors,
    guide = guide_legend(nrow = 23)
  )
#dev.off()



