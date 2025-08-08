## TransBiota Project - Transfeminine microbiota analysis
## Calculate differential abundance of genera and taxa with MaAsLin2
## Elaborated by: Jorge Rojas-Vargas


## Load libraries

# Explore microbiome
library(Maaslin2)


## Read data

load("abund_TF_genus.RData")
load("abund_TF_taxa.RData")
load("metadata.RData")

metadata <- read_tsv("/Users/jorgerv/Documents/transbiota/241019_with_silva138_speciateit2/metadata_final.txt")

# Convert to data frames
abund_TF_genus <- as.data.frame(abund_TF_genus)
abund_TF_taxa <- as.data.frame(abund_TF_taxa)
metadata.maaslin2 <- as.data.frame(metadata)

abund_TF_genus <- abund_TF_genus[apply(abund_TF_genus, 1, function(row) any(row > 0.01)), ]
abund_TF_taxa <- abund_TF_taxa[apply(abund_TF_taxa, 1, function(row) any(row > 0.01)), ]


# Rownames of metadata and subsample_grande
rownames(metadata.maaslin2) <- metadata.maaslin2$UID
metadata.maaslin2$BC <- as.factor(metadata.maaslin2$BC)
metadata.maaslin2$t0_q16 <- as.factor(metadata.maaslin2$t0_q16)
metadata.maaslin2$t0_q16a <- as.factor(metadata.maaslin2$t0_q16a)
metadata.maaslin2$t1_q16 <- as.factor(metadata.maaslin2$t1_q16)


## ----------------------------------------
##
## Run the analysis
##
## ----------------------------------------


## ----------------------------------------
## At genus level with LM models
## Used primarily for continuous data that follows a normal distribution
## ----------------------------------------

## Using abundance tables

# W/linear model
# DEFAULT: normalization=TSS, transform=LOG, method=LM, correction=BH

# Behavioral clusters 
fit_data <- Maaslin2(
  input_data = t(abund_TF_genus),
  input_metadata = metadata.maaslin2,
  output = "maaslind2_TF_t0q16_abund_genus",
  fixed_effects = c("t0_q16"),
  random_effects = c("participant"),
  reference = "t0_q16,0"
)
## **Key Findings from the previous analysis**:
## - There are no associations to plot
## Not even taking as references BC_2, BC_3, or BC_4

# Circumcision status
fit_data <- Maaslin2(
  input_data = t(abund_TF_genus),
  input_metadata = metadata.maaslin2,
  output = "maaslind2_TF_CS_abund_genus",
  fixed_effects = c("CS"),
  random_effects = c("participant")
)
## **Key Findings from the previous analysis**:
## Significant effects with 10 genera!

# Symptoms
fit_data = Maaslin2(
  input_data = t(abund_TF_genus),
  input_metadata = metadata.maaslin2,
  output = "maaslind2_TF_symptoms_abund_genus",
  fixed_effects = c("Symptom"),
  random_effects = c("participant"),
  reference = "Symptom,None")
## **Key Findings from the previous analysis**:
## Significant effects with 20 genera!



## ----------------------------------------
## At taxa level with LM models
## Used primarily for continuous data that follows a normal distribution
## ----------------------------------------

## Using abundance tables

# W/linear model
# DEFAULT: normalization=TSS, transform=LOG, method=LM, correction=BH

# Behavioral clusters 
fit_data <- Maaslin2(
  input_data = t(abund_TF_taxa),
  input_metadata = metadata.maaslin2,
  output = "maaslind2_TF_BC1_abund_taxa",
  fixed_effects = c("BC"),
  random_effects = c("participant"),
  reference = "BC,BC_4"
)
## **Key Findings from the previous analysis**:
## Significant effects with 3 taxa!

# Circumcision status
fit_data <- Maaslin2(
  input_data = t(abund_TF_taxa),
  input_metadata = metadata.maaslin2,
  output = "maaslind2_TF_CS_abund_taxa",
  fixed_effects = c("CS"),
  random_effects = c("participant")
)
## **Key Findings from the previous analysis**:
## Significant effects with 8 taxa!

# Symptoms
fit_data = Maaslin2(
  input_data = t(abund_TF_taxa),
  input_metadata = metadata.maaslin2,
  output = "maaslind2_TF_symptoms_abund_taxa",
  fixed_effects = c("Symptom"),
  random_effects = c("participant"),
  reference = "Symptom,None")
## **Key Findings from the previous analysis**:
## Significant effects with 37 taxa!


