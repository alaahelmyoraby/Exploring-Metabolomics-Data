```markdown
# Exploring Metabolomics Data

This document outlines a workflow for exploring and analyzing metabolomics data, including PCA, normalization, outlier removal, regression analysis, and significant feature identification.

---

## Prerequisites

Ensure that R is installed and the following libraries are loaded:
```r
# Load libraries for visualization and data manipulation
library(ggplot2)
library(dplyr)
```

---

## Steps

### 1. Set Up Working Directory
Set the working directory to the location of your data files:
```r
# Set the directory where data files are located
setwd("C:/My pc/Egcombio/IOSB/lec")
```

### 2. Load Datasets
Load the phenotype and metabolic datasets:
```r
# Load phenotype data
phenotype <- read.csv("Phenotype_synthetic.csv")

# Load metabolic data
metabolic <- read.csv("metabolic.csv")
```

### 3. Check Dimensions of Datasets
Ensure the datasets have been loaded correctly:
```r
# Print dimensions of phenotype and metabolic data
dim(phenotype)
dim(metabolic)
```

### 4. Match and Reorder Datasets
Align phenotype and metabolic data by matching their sample identifiers:
```r
# Match metabolic data samples with phenotype data
ord <- match(metabolic[,"Sample"], phenotype[,"s"])

# Reorder datasets based on matching order
phenotype <- phenotype[ord[which(!is.na(ord))], ]
metabolic <- metabolic[ord[which(!is.na(ord))], ]
```

### 5. Perform PCA on Raw Data
Exclude non-numeric columns and conduct PCA:
```r
# Exclude identifier column for PCA
metab2 <- metabolic[, -1]

# Perform PCA on raw data
pc <- prcomp(metab2)

# Plot PCA eigenvalues
barplot(pc$sdev, main = "PCA Eigenvalues (Raw)")

# Create PCA score plot
plot(pc$x[, 1], pc$x[, 2], main = "PCA Score Plot (Raw)", xlab = "PC1", ylab = "PC2")
```

### 6. Normalize Data Using Z-Scores
Normalize the data to remove scale differences:
```r
# Define Z-score normalization function
Zscore <- function(M) {
  return((M - mean(M, na.rm = TRUE)) / sd(M, na.rm = TRUE))
}

# Apply Z-score normalization to all columns
metab4 <- metab2
for (i in 1:ncol(metab4)) {
  metab4[, i] <- Zscore(metab4[, i])
}
```

### 7. Remove Outliers
Define and apply a function to handle outliers:
```r
# Define function to replace outliers beyond mean ± 3*SD with NA
replace_outliers <- function(x) {
  m <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  lower_bound <- m - 3 * sd_x
  upper_bound <- m + 3 * sd_x
  x[x < lower_bound | x > upper_bound] <- NA
  return(x)
}

# Apply outlier removal to all columns
metab4 <- metab4 %>% mutate(across(everything(), replace_outliers))
```

### 8. Filter Metabolites
Remove metabolites with more than 50% missing values:
```r
# Retain columns with <50% missing values
metab4 <- metab4[, which(colMeans(!is.na(metab4)) > 0.5)]
```

### 9. Perform Log Transformation and Recheck PCA
Log-transform the data to address skewness:
```r
# Log-transform the data
metab_log <- log(metab2 + 1)

# Perform PCA on log-transformed data
pca_log <- prcomp(metab_log, scale. = TRUE)

# Visualize PCA results
ggplot(as.data.frame(pca_log$x), aes(x = PC1, y = PC2)) +
  geom_point() +
  labs(title = "PCA Plot (Log Transformed)", x = "PC1", y = "PC2") +
  theme_dark()
```

### 10. Perform Regression Analysis
Analyze the association between phenotypes and metabolites:
```r
# Initialize a results matrix to store regression results
results_matrix <- matrix(0, nrow = ncol(metab4), ncol = 6)
rownames(results_matrix) <- paste("Metabolite", 1:ncol(metab4))
colnames(results_matrix) <- c("Beta_NoCovariate", "SE_NoCovariate", "Pval_NoCovariate",
                              "Beta_WithCovariate", "SE_WithCovariate", "Pval_WithCovariate")

# Loop through each metabolite and perform regression
for (i in 1:ncol(metab4)) {
  # Linear regression without covariates
  rs1 <- lm(metab4[, i] ~ phenotype[,"Phenotype1"])

  # Linear regression with covariates
  rs2 <- lm(metab4[, i] ~ phenotype[,"Phenotype1"] + phenotype[,"Phenotype1.1"])
  
  # Store results for each regression
  results_matrix[i, ] <- c(
    summary(rs1)$coefficients[2, 1], summary(rs1)$coefficients[2, 2], summary(rs1)$coefficients[2, 4],
    summary(rs2)$coefficients[2, 1], summary(rs2)$coefficients[2, 2], summary(rs2)$coefficients[2, 4]
  )
}
```

### 11. Apply Bonferroni Correction
Identify significant metabolites after multiple testing correction:
```r
# Calculate Bonferroni-corrected threshold
alpha <- 0.05
num_tests <- nrow(results_matrix)
bonferroni_threshold <- alpha / num_tests

# Identify significant metabolites
results_df <- as.data.frame(results_matrix)
significant_metabolites <- results_df[results_df$Pval_NoCovariate < bonferroni_threshold, ]
significant_metabolites <- significant_metabolites[order(significant_metabolites$Pval_NoCovariate), ]
```

### 12. Save Results
Save significant metabolites to a CSV file:
```r
# Save significant metabolites to file
write.csv(significant_metabolites, "significant_metabolites.csv", row.names = TRUE)
```

---

## Outputs
- **PCA Visualizations**: Score plots and eigenvalues for raw, normalized, and log-transformed data.
- **Filtered Dataset**: Processed metabolomics data with outliers and high-missing-value columns removed.
- **Significant Metabolites**: A CSV file (`significant_metabolites.csv`) listing features significantly associated with the phenotype.
