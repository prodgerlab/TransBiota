## Load libraries

# Statistical tools
library(rstatix)

# Handle tables
library(dplyr)
library(tibble)

# Plot tools
library(ggplot2)
library(reshape2)
library(svglite)

#load in xlxs files

#make first column of abund_TF_genus header
abund_TF_genus <- column_to_rownames(abund_TF_genus, var = "...1")

abund_TF_genus_filtered <- abund_TF_genus[rownames(abund_TF_genus) != "Unassigned", , drop = FALSE] 

abund_TF_t <- as.data.frame(t(abund_TF_genus_filtered))
rownames(metadata_nugent) <- metadata_nugent$UID
combined_data <- cbind(metadata_nugent, abund_TF_t)

# Define the bacteria c("T1", "T2", "t3")
TF <- c("Rikenellaceae_RC9_gut_group", "Staphylococcus", "Kallipyga", "Anaeroglobus", "Aerococcus", "Peptostreptococcus", "Peptococcus", "Negativicoccus", "W5053", "Facklamia", "S5-A14a", "Parvimonas", "Fusobacterium", "Streptococcus", "Ezakiella", "Finegoldia", "Anaerococcus", "Peptoniphilus", "Campylobacter", "Varibaculum", "Mobiluncus", "Gleimia", "Mogibacterium", "Slackia", "Moryella", "Actinotignum", "Howardella", "Corynebacterium", "Propionimicrobium", "Atopobium", "Fastidiosipila", "Arcanobacterium", "Murdochiella", "Schaalia", "Fenollaria", "Prevotella", "Dialister", "Porphyromonas", "Hoylesella", "Fannyhessea", "Gardnerella", "Bifidobacterium", "Lawsonella", "Lactobacillus")

# Define the columns of interest
cytokine_cols <- c("IL1A", "IL1B", "IL6", "IL8", "MIG", "MIP1B", "RANTES")
columns_of_interest <- c(cytokine_cols, TF)

combined_data_selected <- combined_data[, columns_of_interest]

# Compute Spearman's correlation
set.seed(123)
cormat <- combined_data_selected %>% cor_mat(method = "spearman", alternative = "two.sided", conf.level = 0.95)

# Get p-values for Spearman's correlation
cormat.p <- cormat %>% cor_get_pval()

cormat <- column_to_rownames(cormat, var = "rowname")
cormat.p <- column_to_rownames(cormat.p, var = "rowname")

# Taxa vs Cytokines Spearman's correlation
cormat.imm <- cormat[TF, cytokine_cols]
cormat.p.imm <- cormat.p[TF, cytokine_cols]

# Convert to matrix
cormat <- as.matrix(cormat)
cormat.imm <- as.matrix(cormat.imm)

# Melt the correlation matrix
melted_cormat <- melt(cormat, na.rm = TRUE)
melted_cormat.imm <- melt(cormat.imm, na.rm = TRUE)

# Round the values to 2 decimal places
melted_cormat <- melted_cormat %>%
  dplyr::mutate(value = round(value, 2))
melted_cormat.imm <- melted_cormat.imm %>% 
  dplyr::mutate(value = round(value, 2))

# Rename the columns in both the immune subset dataset and the original bacteria correlation dataset 
# so they can be matched
names(melted_cormat.imm)[names(melted_cormat.imm) == "Var1"] <- "taxa"
names(melted_cormat)[names(melted_cormat) == "Var2"] <- "taxa"

# Convert the taxa column to a factor with desired order
melted_cormat$taxa <- factor(melted_cormat$taxa, levels = TF)
melted_cormat.imm$taxa <- factor(melted_cormat.imm$taxa, levels = TF)

# Use the inner_join function to match the two datasets according to the taxa name order
order <- melted_cormat %>%
  dplyr::select(taxa) %>%
  inner_join(melted_cormat.imm, by = "taxa")

# Ensure there are no duplicates
order <- distinct(order)

# Create immune-to-bacteria correlation heatmap
heatmap_immune_strip <- ggplot(order, aes(x = Var2, y = taxa, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#01baef", high = "#8ac926", mid = "white", 
                       midpoint = 0, limit = c(-0.6, 0.6), space = "Lab",
                       name="Spearman\nCorrelation") +
  theme_minimal() + # Minimal theme
  theme(
    axis.text.y = element_text(size = 7, vjust= 1, hjust = 1),
    axis.text.x = element_text(size = 7, angle = 45, vjust= 1, hjust = 1)) +
  scale_x_discrete(labels = cytokine_cols) +
  coord_fixed(ratio=0.75) + 
  geom_text(aes(Var2, taxa, label = value), color = "black", size = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())


heatmap_immune_strip

# p values
print(cormat.p.imm)


#To make Prevalence heatmap 


library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(scales)


abund_TF_genus_filtered <- abund_TF_genus[TF,]

abund_TF_genus_filtered <- abund_TF_genus_filtered %>% 
  tibble::rownames_to_column(var = "Bacteria")

print(abund_TF_genus_filtered)

# Convert to long format
data_long <- abund_TF_genus_filtered %>%
  pivot_longer(cols = -1, names_to = "Participant", values_to = "Relative_Abundance")
head(data_long)

# Compute prevalence: proportion of participants where the bacterium is present (>0)
prevalence_data <- data_long %>%
  group_by(Bacteria) %>%
  summarise(Prevalence = sum(Relative_Abundance > 0, na.rm = TRUE) / n())

# Check prevalence table
head(prevalence_data)

# Define the custom order inside ggplot using factor()
custom_order <- c("Rikenellaceae_RC9_gut_group", "Staphylococcus", "Kallipyga", "Anaeroglobus", "Aerococcus", "Peptostreptococcus", "Peptococcus", "Negativicoccus", "W5053", "Facklamia", "S5-A14a", "Parvimonas", "Fusobacterium", "Streptococcus", "Ezakiella", "Finegoldia", "Anaerococcus", "Peptoniphilus", "Campylobacter", "Varibaculum", "Mobiluncus", "Gleimia", "Mogibacterium", "Slackia", "Moryella", "Actinotignum", "Howardella", "Corynebacterium", "Propionimicrobium", "Atopobium", "Fastidiosipila", "Arcanobacterium", "Murdochiella", "Schaalia", "Fenollaria", "Prevotella", "Dialister", "Porphyromonas", "Hoylesella", "Fannyhessea", "Gardnerella", "Bifidobacterium", "Lawsonella", "Lactobacillus")  # Replace with actual bacteria names

#make graph
prevalence_heatmap <- ggplot(prevalence_data, aes(x = "", y = factor(Bacteria, levels = custom_order), fill = Prevalence)) +
  geom_tile(color = "black") +  # Add borders around tiles
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1), name = "Prevalence") +
  geom_text(aes(label = scales::percent(Prevalence, accuracy = 1)), color = "black", size = 1.75) +  # Add prevalence number inside tiles
  theme_minimal() +
  labs(
    title = "Bacterial Prevalence Heatmap",
    x = "",
    y = "Bacteria"
  ) +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis text 
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10) 
  ) +
  coord_fixed(ratio = 1) 

print(prevalence_heatmap)

print(prevalence_data, n=44)



#To make relative abundance barplot

library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(scales)


abund_TF_genus_filtered <- abund_TF_genus[TF,]

abund_TF_genus_filtered <- abund_TF_genus_filtered %>% 
  tibble::rownames_to_column(var = "Bacteria")

# Convert to long format
data_long <- abund_TF_genus_filtered %>%
  pivot_longer(cols = -1, names_to = "Participant", values_to = "Relative_Abundance")
head(data_long)

# Compute median and quartiles excluding zero values
median_data <- data_long %>%
  filter(Relative_Abundance > 0) %>%  # Exclude zero values
  group_by(Bacteria) %>%
  summarise(
    Median_Abundance = median(Relative_Abundance, na.rm = TRUE),
    Q1 = quantile(Relative_Abundance, 0.25, na.rm = TRUE),  # First quartile
    Q3 = quantile(Relative_Abundance, 0.75, na.rm = TRUE)   # Third quartile
  )

# Define the custom order inside ggplot using factor()
custom_order <- c("Rikenellaceae_RC9_gut_group", "Staphylococcus", "Kallipyga", "Anaeroglobus", "Aerococcus", "Peptostreptococcus", "Peptococcus", "Negativicoccus", "W5053", "Facklamia", "S5-A14a", "Parvimonas", "Fusobacterium", "Streptococcus", "Ezakiella", "Finegoldia", "Anaerococcus", "Peptoniphilus", "Campylobacter", "Varibaculum", "Mobiluncus", "Gleimia", "Mogibacterium", "Slackia", "Moryella", "Actinotignum", "Howardella", "Corynebacterium", "Propionimicrobium", "Atopobium", "Fastidiosipila", "Arcanobacterium", "Murdochiella", "Schaalia", "Fenollaria", "Prevotella", "Dialister", "Porphyromonas", "Hoylesella", "Fannyhessea", "Gardnerella", "Bifidobacterium", "Lawsonella", "Lactobacillus")  # Replace with actual bacteria names

# Bar plot with Q1-Q3 error bars
ggplot(median_data, aes(x = factor(Bacteria, levels = custom_order), y = Median_Abundance)) +
  geom_bar(stat = "identity", fill = "snow2", color = "black") +  # Create bars
  geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.3) +  # Error bars using Q1 & Q3
  theme_minimal() +
  labs(
    title = "Median Relative Abundance of Bacteria (Excluding Zeros)",
    x = "Bacteria",
    y = "Median Relative Abundance"
  ) +
  theme(
    axis.text.x = element_text(hjust = 1),  # Rotate x-axis labels
    panel.grid.major = element_line(color = "white"),  # Make major grid lines white
    panel.grid.minor = element_line(color = "white")   # Make minor grid lines white
  ) +
  coord_flip()

print(median_data, n=44)



#To make Nugent by bacteria heatmap

# Statistical tools
library(rstatix)

# Handle tables
library(dplyr)
library(tibble)

# Plot tools
library(ggplot2)
library(reshape2)
library(svglite)

# Filter out "Unassigned" taxa
abund_TF_genus_filtered <- abund_TF_genus[rownames(abund_TF_genus) != "Unassigned", , drop = FALSE]

# Transpose taxa data and combine with metadata_nugent
abund_TF_t <- as.data.frame(t(abund_TF_genus_filtered))
rownames(metadata_nugent) <- metadata_nugent$UID
combined_data <- cbind(metadata_nugent, abund_TF_t)

# Define taxa of interest
TF <- c("Rikenellaceae_RC9_gut_group", "Staphylococcus", "Kallipyga", "Anaeroglobus", 
        "Aerococcus", "Peptostreptococcus", "Peptococcus", "Negativicoccus", "W5053", 
        "Facklamia", "S5-A14a", "Parvimonas", "Fusobacterium", "Streptococcus", 
        "Ezakiella", "Finegoldia", "Anaerococcus", "Peptoniphilus", "Campylobacter", 
        "Varibaculum", "Mobiluncus", "Gleimia", "Mogibacterium", "Slackia", 
        "Moryella", "Actinotignum", "Howardella", "Corynebacterium", "Propionimicrobium", 
        "Atopobium", "Fastidiosipila", "Arcanobacterium", "Murdochiella", "Schaalia", 
        "Fenollaria", "Prevotella", "Dialister", "Porphyromonas", "Hoylesella", 
        "Fannyhessea", "Gardnerella", "Bifidobacterium", "Lawsonella", "Lactobacillus")

# Select columns of interest: Nugent score + taxa
nugent_col <- "Nugent"
columns_of_interest <- c(nugent_col, TF)
combined_data_selected <- combined_data[, columns_of_interest]

# Calculate Spearman correlation matrix and p-values
set.seed(123)
cormat <- combined_data_selected %>% cor_mat(method = "spearman", alternative = "two.sided", conf.level = 0.95)
cormat.p <- cormat %>% cor_get_pval()

# Convert to data.frames with rownames as columns for easy subsetting
cormat <- column_to_rownames(cormat, var = "rowname")
cormat.p <- column_to_rownames(cormat.p, var = "rowname")

# Extract bacteria vs Nugent correlations and p-values
cormat.nug <- cormat[TF, nugent_col, drop = FALSE]
cormat.p.nug <- cormat.p[TF, nugent_col, drop = FALSE]

# Melt the bacteria vs Nugent correlation matrix for plotting
melted_cormat.nug <- melt(as.matrix(cormat.nug), na.rm = TRUE)

# Fix factor levels and names to avoid NAs
melted_cormat.nug <- melted_cormat.nug %>% 
  rename(taxa = Var1, Var2 = Var2) %>% 
  mutate(
    taxa = as.character(taxa),
    taxa = factor(taxa, levels = TF),       
    Var2 = factor(Var2, levels = "Nugent"),    
    value = round(value, 2)
  )

# Should show taxa names, Var2 = Nugent, and values
print(head(melted_cormat.nug))

# Plot heatmap
heatmap_nugent_strip <- ggplot(melted_cormat.nug, aes(x = Var2, y = taxa, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = value), color = "black", size = 2) +
  scale_fill_gradient2(
    low = "#A8D0FF", mid = "white", high = "#FFC0CB",
    midpoint = 0, limit = c(-0.5, 0.5), space = "Lab",
    name = "Spearman\nCorrelation"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed(ratio = 0.75) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank()
  )

print(heatmap_nugent_strip)

# print p-values for inspection
print(cormat.p.nug)

