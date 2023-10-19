###Spearman with p-signifiance

# Load required libraries
library(dplyr)
library(igraph)  

setwd("/Users/user/Documents/Metagenomics_paper/Correlation")

# Read the genus and KEGG pathway data
genus_data <- read.csv("Genus.csv", header = TRUE, row.names = 1)
kegg_data <- read.csv("KEGG.csv", header = TRUE, row.names = 1)

# Normalize the data (convert to proportions)
genus_data_normalized <- genus_data / rowSums(genus_data)
kegg_data_normalized <- kegg_data / rowSums(kegg_data)

# Convert normalized data frames to numeric matrices
genus_data_numeric <- as.matrix(genus_data_normalized)
kegg_data_numeric <- as.matrix(kegg_data_normalized)

# Get the number of columns in the data matrices
num_cols_genus <- ncol(genus_data_numeric)
num_cols_kegg <- ncol(kegg_data_numeric)

# Initialize matrices to store the correlation coefficients and p-values
cor_matrix <- matrix(NA, nrow = num_cols_genus, ncol = num_cols_kegg)
p_values <- matrix(NA, nrow = num_cols_genus, ncol = num_cols_kegg)

# Calculate Spearman correlations and p-values for all pairs of columns
for (i in 1:num_cols_genus) {
  for (j in 1:num_cols_kegg) {
    cor_test_result <- cor.test(genus_data_numeric[, i], kegg_data_numeric[, j], method = "spearman")
    cor_matrix[i, j] <- cor_test_result$estimate
    p_values[i, j] <- cor_test_result$p.value
  }
}

# Apply a significance threshold (e.g., p < 0.05)
significant_cor_matrix <- cor_matrix * (p_values < 0.05)

# Export the significant correlation matrix as a CSV file
write.csv(significant_cor_matrix, "genus_kegg_significant_correlations.csv")