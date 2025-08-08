## TransBiota Project - Transfeminine microbiota analysis
## Obtain Spearman correlation heatmap
## Elaborated by: Jorge Rojas-Vargas

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

## Read data
load("abund_TF_genus.RData")
load("abund_CM_genus.RData")
load("metadata.RData")

## --------------------------------------------------------------------------------
##
## Spearman correlation and Taxa Clusters definition for TF cohort
##
## --------------------------------------------------------------------------------

abund_TF_genus_filtered <- abund_TF_genus[rownames(abund_TF_genus) != "Unassigned", , drop = FALSE]
medians_TF <- apply(abund_TF_genus_filtered, 1, median)
top30_TF <- names(sort(medians_TF, decreasing = TRUE)[1:30])

abund_TF_t <- as.data.frame(t(abund_TF_genus))
rownames(metadata) <- metadata$UID

samples_in_both <- intersect(rownames(metadata), rownames(abund_TF_t))

metadata <- metadata[samples_in_both, , drop = FALSE]
abund_TF_t <- abund_TF_t[samples_in_both, , drop = FALSE]

combined_data <- cbind(metadata, abund_TF_t)

# Columns of interest
cytokine_cols <- c("logIL1A", "logIL1B", "logIL6", "logIL8", "logMIG", "logMIP1B", "logRANTES")
#cytokine_cols <- c("IL1A", "IL1B", "IL6", "IL8", "MIG", "MIP1B", "RANTES")
columns_of_interest <- c(cytokine_cols, top30_TF)

# Select columns
combined_data_selected <- combined_data[, columns_of_interest]

# Calculate Spearman's correlation
set.seed(123)
cormat <- combined_data_selected %>% cor_mat(method = "spearman", alternative = "two.sided", conf.level = 0.95)

# Get p-values for Spearman's correlation
cormat.p <- cormat %>% cor_get_pval()

cormat <- column_to_rownames(cormat, var = "rowname")
cormat.p <- column_to_rownames(cormat.p, var = "rowname")

# Taxa vs Cytokines Spearman's correlation
cormat.imm <- cormat[top30_TF, cytokine_cols]
cormat.p.imm <- cormat.p[top30_TF, cytokine_cols]

# Reorder rows and columns of cormat based on the order of top_TF_genus
cormat <- cormat[top30_TF, top30_TF]
cormat.p <- cormat.p[top30_TF, top30_TF]

# Benjamini-Hochberg (BH) or False Discovery Rate adjustment
# Adjust p-values using the BH method for subsamples
fdrs_cytokines <- p.adjust(as.matrix(cormat.p.imm), method = "BH")
fdrs_bacteria <- p.adjust(as.matrix(cormat.p), method = "BH")
# Convert the adjusted p-values back to the matrix structure
fdrs_matrix_cytokines <- matrix(fdrs_cytokines, nrow = nrow(cormat.p.imm), ncol = ncol(cormat.p.imm))
rownames(fdrs_matrix_cytokines) <- rownames(cormat.p.imm)
colnames(fdrs_matrix_cytokines) <- colnames(cormat.p.imm)
fdrs_matrix_bacteria <- matrix(fdrs_bacteria, nrow = nrow(cormat.p), ncol = ncol(cormat.p))
rownames(fdrs_matrix_bacteria) <- rownames(cormat.p)
colnames(fdrs_matrix_bacteria) <- colnames(cormat.p)

## Clustering with k-means

# Define a function to calculate the SSE for different k values
calculate_sse <- function(data, max_k = 10) {
  sse <- numeric(max_k)
  for (k in 1:max_k) {
    set.seed(123)
    kmeans_result <- kmeans(data, centers = k)
    sse[k] <- sum(kmeans_result$withinss)
  }
  return(sse)
}

# Calculate the SSE for k values from 1 to 10
max_k <- 20
sse <- calculate_sse(cormat, max_k)
sse
plot(1:max_k, sse, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters (k)",
     ylab = "Sum of Squared Errors (SSE)",
     main = "Elbow Method for Optimal k")

# Perform clustering with k-means
k <- 4  # Initial number of clusters (this number can be adjusted)
set.seed(123)
km_result <- kmeans(cormat, centers = k)
km_result$cluster

# Ensure each cluster has at least 3 taxa
min_taxa_per_cluster <- 3

# Function to ensure the cluster size constraint
adjust_clusters <- function(km_result, min_taxa_per_cluster, cormat) {
  cluster_sizes <- table(km_result$cluster)
  
  # Find clusters that do not meet the constraint
  small_clusters <- which(cluster_sizes < min_taxa_per_cluster)
  
  while(length(small_clusters) > 0) {
    for (small_cluster in small_clusters) {
      # Get taxa from the small cluster
      taxa_in_small_cluster <- which(km_result$cluster == small_cluster)
      
      for (taxa in taxa_in_small_cluster) {
        # Reassign taxa to the closest cluster that meets the constraint
        dist_to_other_clusters <- sapply(setdiff(1:k, small_cluster), function(cl) {
          mean(as.dist(cormat[taxa, km_result$cluster == cl]))
        })
        
        new_cluster <- which.min(dist_to_other_clusters)
        km_result$cluster[taxa] <- new_cluster
      }
    }
    
    cluster_sizes <- table(km_result$cluster)
    small_clusters <- which(cluster_sizes < min_taxa_per_cluster)
  }
  
  return(km_result)
}

# Adjust clusters to meet the constraint
adjusted_clusters <- adjust_clusters(km_result, min_taxa_per_cluster, cormat)
taxa_clusters <- adjusted_clusters$cluster
taxa_clusters

# Reorder the correlation matrix based on the adjusted clusters
reorder_cormat_by_clusters <- function(cormat, clusters) {
  order_idx <- order(clusters)
  cormat <- cormat[order_idx, order_idx]
  return(cormat)
}

cormat_reordered <- reorder_cormat_by_clusters(cormat, adjusted_clusters$cluster)

cormat <- cormat_reordered

# Plot
cormat <- as.matrix(cormat)
melted_cormat <- melt(cormat, na.rm = TRUE)

# Round the values to 2 decimal places
melted_cormat <- melted_cormat %>%
  dplyr::mutate(value = round(value, 2))

ggheatmap <- ggplot(melted_cormat, aes(x = factor(Var2, levels = rev(levels(Var2))), y = Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#383da0", high = "#d32b07", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 11, hjust = 1), 
        axis.text.y = element_text(size = 15, hjust = 1))+
  labs(x=NULL, y= NULL) +
  coord_fixed() + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 3.5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

ggheatmap

cormat.imm <- as.matrix(cormat.imm)

# Melt the correlation matrix
melted_cormat.imm <- melt(cormat.imm, na.rm = TRUE)

# Reduce Spearman correlation values to 2 sig digs
melted_cormat.imm<-melted_cormat.imm %>% 
  dplyr::mutate(value = round(value, 2))

# Rename the columns in both the immune subset data set and the original bacteria correlation data set so I can match them up 
names(melted_cormat.imm)[names(melted_cormat.imm) == "Var1"] <- "taxa"
names(melted_cormat)[names(melted_cormat) == "Var2"] <- "taxa"

# Use the inner_join function for match the two data sets up according to the taxa name order that was determined by the reordering function above 
order<-melted_cormat%>%
  dplyr::select(taxa)%>%
  inner_join(melted_cormat.imm, by = "taxa")

# Ensure there are no duplicates
order <- distinct(order)

# The 'order' data set is ready to be plotted

# Create immune to bacteria correlation heat map
heatmap_immune_strip<-ggplot(order, aes(x = Var2, y = taxa, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#01baef", high = "#8ac926", mid = "white", 
                       midpoint = 0, limit = c(-0.5, 0.5), space = "Lab",
                       name="Spearman\nCorrelation") +
  theme_minimal() + # minimal theme
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11, angle = 45, vjust= 1, hjust = 1)) +
  scale_x_discrete(labels = cytokine_cols) +
  coord_fixed() + 
  geom_text(aes(Var2, taxa, label = value), color = "black", size = 3.5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

heatmap_immune_strip

svglite("Fig_4A.svg", width = 22, height = 17)
(ggheatmap | heatmap_immune_strip)
dev.off()





# Assign clusters based on genus name
cluster_map <- list(
  TC1 = "Dialister",
  TC2 = "Prevotella",
  TC3 = "Streptococcus",
  TC4 = "Hoylesella"
)

# Filter the genera corresponding to each TC
TC1_taxa <- names(taxa_clusters)[taxa_clusters == taxa_clusters["Dialister"]]
TC2_taxa <- names(taxa_clusters)[taxa_clusters == taxa_clusters["Prevotella"]]
TC3_taxa <- names(taxa_clusters)[taxa_clusters == taxa_clusters["Streptococcus"]]
TC4_taxa <- names(taxa_clusters)[taxa_clusters == taxa_clusters["Hoylesella"]]

# Print the taxa cluster members
TC1_taxa
TC2_taxa
TC3_taxa
TC4_taxa

# Create the abund_TF_cluster dataframe with columns for TC1, TC2, TC3, TC4, and "Other"
abund_TF_cluster <- data.frame(
  TC1 = rowSums(abund_TF_t[, TC1_taxa, drop = FALSE]),
  TC2 = rowSums(abund_TF_t[, TC2_taxa, drop = FALSE]),
  TC3 = rowSums(abund_TF_t[, TC3_taxa, drop = FALSE]),
  TC4 = rowSums(abund_TF_t[, TC4_taxa, drop = FALSE]),
  Other = rowSums(abund_TF_t[, !(colnames(abund_TF_t) %in% top30_TF), drop = FALSE])
)

# Verify the result
head(abund_TF_cluster)


## Combining taxa clusters with metadata

mdt_table <- as.data.frame(metadata)
rownames(mdt_table) <- metadata$UID
mdt_table$UID <- NULL

samples_in_both <- intersect(rownames(mdt_table), rownames(abund_TF_cluster))

mdt_table_filtered <- mdt_table[samples_in_both, , drop = FALSE]
abund_TF_cluster_filtered <- abund_TF_cluster[samples_in_both, , drop = FALSE]

# If using the TCs abundance only
combined_data_cluster <- cbind(mdt_table_filtered, abund_TF_cluster_filtered)




## --------------------------------------------------------------
##
## Stackbar plots to get sample clusters (Fig. 4A)
##
## --------------------------------------------------------------

##### Stack bar plots for TC to find the Sample Cluster

# Subsample with abundances
taxa_data <- abund_TF_cluster
taxa_data$UID <- rownames(taxa_data)

# Organize according to hierarchical clustering

# Calculate distances between samples
dist_mat <- dist(taxa_data[, colnames(taxa_data)[1:5]], method = "euclidean")

# Apply hierarchical clustering to order the samples
clusters <- hclust(dist_mat, method = "ward.D2")


# Make cut with cutree
h_value = 1.1 # 6 clusters
clus <- cutree(clusters, h = h_value)

# Add the cluster information to the original dataframe
taxa_data$Cluster <- factor(clus)
taxa_data$Cluster_letter <- LETTERS[as.numeric(taxa_data$Cluster)]
cluster_summary <- taxa_data %>%
  group_by(Cluster_letter) %>%
  summarise(Number_members = n())
taxa_data <- merge(taxa_data, cluster_summary, by = "Cluster_letter")
taxa_data$Final_cluster <- paste(taxa_data$Number_members, taxa_data$Cluster_letter, sep = "")

# Separate the numeric and letter components of `Final_cluster`
taxa_data_ordered <- taxa_data %>%
  mutate(
    Number = as.numeric(gsub("[A-Z]", "", Final_cluster)),  # Extract the number
    Letter = gsub("[0-9]", "", Final_cluster)  # Extract the letter
  ) %>%
  arrange(desc(Number), Letter, UID)  # Order by number (descending), then by letter, and finally by UID

# Create a df from taxa_data_ordered
sample_TGC <- taxa_data_ordered[, c("UID", "Final_cluster")]

# Create a factor based on the unique values of Final_cluster, ordered in descending order
sample_TGC$Final_cluster_factor <- factor(sample_TGC$Final_cluster, levels = unique(sample_TGC$Final_cluster[order(-as.numeric(gsub("[A-Z]", "", sample_TGC$Final_cluster)))]))

# Create a TGC column where letters are assigned in ascending order but starting from the highest value
sample_TGC$TGC <- LETTERS[as.numeric(sample_TGC$Final_cluster_factor)]

combined_data_cluster_w_TGC <- combined_data_cluster
combined_data_cluster_w_TGC$UID <- rownames(combined_data_cluster)
combined_data_cluster_w_TGC <- merge(combined_data_cluster_w_TGC, sample_TGC[, c("UID", "TGC")], by = "UID", all.x = TRUE)

# Remove auxiliary columns
taxa_data_ordered$Cluster <- NULL
taxa_data_ordered$Cluster_letter <- NULL
taxa_data_ordered$Number_members <- NULL
taxa_data_ordered$Final_cluster <- NULL
taxa_data_ordered$Letter <- NULL
taxa_data_ordered$Number <- NULL

## -----------------
## Figure 4A
## -----------------

# Convert the data to long format for ggplot2
taxa_data_long <- melt(taxa_data_ordered, id.vars = "UID")

# Define custom colors
custom_colors <- c("#7209b7", "#e63946", "#ff9f1c", "#023e8a","#adb5bd")

# Create the stacked bar plot
#svglite("samples_taxa_cluster_24.11.17.svg", width = 22, height = 10)
ggplot(taxa_data_long, aes(x = factor(UID, levels = taxa_data_ordered$UID), y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Samples", y = "Relative Abundance", fill = "Taxa") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 20),
        legend.position = "right",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 22)) +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  ggtitle("Stacked Bar Plot of Taxa Relative Abundances")
#dev.off()

#svglite("samples_taxa_cluster_dendogram.svg", width = 10, height = 5)
plot(clusters, main = "Dendrogram", hang = -1)
#dev.off()


## ------------------------------------------------------------
##
## Heatmap at genus level in each cluster of samples (Fig. 4B)
##
## ------------------------------------------------------------


# Filter the abund dataframe to include only rows in top_40
abund_top30 <- abund_TF_genus[rownames(abund_TF_genus) %in% top30_TF, ]

# Reorder the rows of abund_top30 according to the order of top_30
abund_top30 <- abund_top30[match(top30_TF, rownames(abund_top30)), ]

# Convert rownames into a column so they can be used as id.vars
abund_top30$Taxa <- rownames(abund_top30)

# Use melt to transform the dataframe into long format
abund_top30_long <- melt(abund_top30, id.vars = "Taxa")

# Vector of desired order for taxa
desired_order <- c(
  "Dialister", "Mobiluncus", "Porphyromonas", "Peptococcus", "W5053", 
  "Fenollaria", "Atopobium", "Prevotella", "Lawsonella", "Parvimonas", "Fusobacterium", 
  "Howardella", "Streptococcus", "Anaerococcus", "Corynebacterium", 
  "Finegoldia", "Schaalia", "Lactobacillus", "Peptoniphilus", "Arcanobacterium", 
  "Propionimicrobium", "Negativicoccus", "Ezakiella", "Actinotignum", 
  "Campylobacter", "Murdochiella", "Hoylesella", "Varibaculum",  
  "S5-A14a", "Fastidiosipila"
)

# Reverse the order of Taxa not including Other
abund_top30_long$Taxa <- factor(abund_top30_long$Taxa, levels = rev(desired_order))
abund_top30_long$variable <- factor(abund_top30_long$variable, levels = taxa_data_ordered$UID)

# Set the range of values for the color scale
value_min <- min(abund_top30_long$value, na.rm = TRUE)
value_max <- max(abund_top30_long$value, na.rm = TRUE)

# Save the heatmap as SVG
#svglite("Fig_4B.svg", width = 22, height = 10)
ggplot(abund_top30_long, aes(x = variable, y = Taxa, fill = value)) +
  geom_tile(color = "white") + # Tiles for heatmap
  theme_minimal() +
  labs(x = "Samples", y = "Taxa", fill = "Relative Abundance") +
  scale_fill_gradientn(
    colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),
    values = scales::rescale(c(value_min, value_max)),
    limits = c(value_min, value_max),
    na.value = "grey50"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 18, hjust = 0),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18)
  ) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 20))
#dev.off()







## --------------------------------------------------------------------------------
##
## Spearman correlation and Taxa Clusters definition for CM cohort
##
## --------------------------------------------------------------------------------

abund_PC_genus_filtered <- abund_PC_genus[rownames(abund_PC_genus) != "Unassigned", , drop = FALSE]
medians_PC <- apply(abund_PC_genus_filtered, 1, median)
top30_PC <- names(sort(medians_PC, decreasing = TRUE)[1:30])

abund_PC_t <- as.data.frame(t(abund_PC_genus))


# Calculate Spearman's correlation
set.seed(123)
cormat <- abund_PC_t %>% cor_mat(method = "spearman", alternative = "two.sided", conf.level = 0.95)

# Get p-values for Spearman's correlation
cormat.p <- cormat %>% cor_get_pval()

cormat <- column_to_rownames(cormat, var = "rowname")
cormat.p <- column_to_rownames(cormat.p, var = "rowname")

# Reorder rows and columns of cormat based on the order of top_PC_genus
cormat <- cormat[top30_PC, top30_PC]
cormat.p <- cormat.p[top30_PC, top30_PC]
write.csv(cormat.p, "Supp_Table_pvalues_CM.csv")

# Benjamini-Hochberg (BH) or False Discovery Rate adjustment
# Adjust p-values using the BH method for subsamples
fdrs_bacteria <- p.adjust(as.matrix(cormat.p), method = "BH")
# Convert the adjusted p-values back to the matrix structure
fdrs_matrix_bacteria <- matrix(fdrs_bacteria, nrow = nrow(cormat.p), ncol = ncol(cormat.p))
rownames(fdrs_matrix_bacteria) <- rownames(cormat.p)
colnames(fdrs_matrix_bacteria) <- colnames(cormat.p)

## Clustering with k-means

# Define a function to calculate the SSE for different k values
calculate_sse <- function(data, max_k = 10) {
  sse <- numeric(max_k)
  for (k in 1:max_k) {
    set.seed(123)
    kmeans_result <- kmeans(data, centers = k)
    sse[k] <- sum(kmeans_result$withinss)
  }
  return(sse)
}

# Calculate the SSE for k values from 1 to 10
max_k <- 20
sse <- calculate_sse(cormat, max_k)
sse
plot(1:max_k, sse, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters (k)",
     ylab = "Sum of Squared Errors (SSE)",
     main = "Elbow Method for Optimal k")

# Perform clustering with k-means
k <- 5  # Initial number of clusters (this number can be adjusted)
set.seed(123)
km_result <- kmeans(cormat, centers = k)
km_result$cluster

# Ensure each cluster has at least 3 taxa
min_taxa_per_cluster <- 3

# Function to ensure the cluster size constraint
adjust_clusters <- function(km_result, min_taxa_per_cluster, cormat) {
  cluster_sizes <- table(km_result$cluster)
  
  # Find clusters that do not meet the constraint
  small_clusters <- which(cluster_sizes < min_taxa_per_cluster)
  
  while(length(small_clusters) > 0) {
    for (small_cluster in small_clusters) {
      # Get taxa from the small cluster
      taxa_in_small_cluster <- which(km_result$cluster == small_cluster)
      
      for (taxa in taxa_in_small_cluster) {
        # Reassign taxa to the closest cluster that meets the constraint
        dist_to_other_clusters <- sapply(setdiff(1:k, small_cluster), function(cl) {
          mean(as.dist(cormat[taxa, km_result$cluster == cl]))
        })
        
        new_cluster <- which.min(dist_to_other_clusters)
        km_result$cluster[taxa] <- new_cluster
      }
    }
    
    cluster_sizes <- table(km_result$cluster)
    small_clusters <- which(cluster_sizes < min_taxa_per_cluster)
  }
  
  return(km_result)
}

# Adjust clusters to meet the constraint
adjusted_clusters <- adjust_clusters(km_result, min_taxa_per_cluster, cormat)
taxa_clusters <- adjusted_clusters$cluster
taxa_clusters

# Reorder the correlation matrix based on the adjusted clusters
reorder_cormat_by_clusters <- function(cormat, clusters) {
  order_idx <- order(clusters)
  cormat <- cormat[order_idx, order_idx]
  return(cormat)
}

cormat_reordered <- reorder_cormat_by_clusters(cormat, adjusted_clusters$cluster)

cormat <- cormat_reordered

# Plot
cormat <- as.matrix(cormat)
melted_cormat <- melt(cormat, na.rm = TRUE)

# Round the values to 2 decimal places
melted_cormat <- melted_cormat %>%
  dplyr::mutate(value = round(value, 2))

ggheatmap <- ggplot(melted_cormat, aes(x = factor(Var2, levels = rev(levels(Var2))), y = Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#383da0", high = "#d32b07", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 11, hjust = 1), 
        axis.text.y = element_text(size = 15, hjust = 1))+
  labs(x=NULL, y= NULL) +
  coord_fixed() + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 3.5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

ggheatmap

svglite("Supp_Fig_8.svg", width = 16, height = 16)
ggheatmap
dev.off()



