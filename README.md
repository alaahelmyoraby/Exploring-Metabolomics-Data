# Exploring-Metabolomics-Data

```markdown
# PCA and Regression Analysis Workflow

This document outlines the steps for performing PCA, data normalization, outlier removal, and regression analysis on phenotype and metabolic datasets.

---

## Prerequisites

- Install R and the following packages:
  - `ggplot2`
  - `dplyr`

## Steps

### 1. Load Necessary Libraries
```r
library(ggplot2)
library(dplyr)
```

### 2. Set Working Directory
Set the working directory to the location of your data files:
```r
setwd("C:/Egcombio/IOSB")
```

### 3. Load Datasets
Load the phenotype and metabolic datasets:
```r
phenotype <- read.csv("Phenotype_synthetic.csv")  # Phenotype data
metabolic <- read.csv("metabolic.csv")  # Metabolic data
```

### 4. Match and Reorder Datasets
Ensure that both datasets are correctly aligned:
```r
ord <- match(metabolic[,"Sample"], phenotype[,"s"])
phenotype <- phenotype[ord[which(!is.na(ord))], ]
metabolic <- metabolic[ord[which(!is.na(ord))], ]
```

### 5. Perform PCA on Raw Data
Exclude non-numeric columns and conduct PCA:
```r
metab2 <- metabolic[, -1]
pc <- prcomp(metab2)
barplot(pc$sdev, main = "PCA Eigenvalues (Raw)")
plot(pc$x[, 1], pc$x[, 2], main = "PCA Score Plot (Raw)", xlab = "PC1", ylab = "PC2")
```

### 6. Normalize Data Using Z-Scores
Apply Z-score normalization to the dataset:
```r
Zscore <- function(M) {
  return((M - mean(M, na.rm = TRUE)) / sd(M, na.rm = TRUE))
}

metab4 <- metab2
for (i in 1:ncol(metab4)) {
  metab4[, i] <- Zscore(metab4[, i])
}
```

### 7. Remove Outliers
Define and apply a function to remove outliers:
```r
replace_outliers <- function(x) {
  m <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  lower_bound <- m - 3 * sd_x
  upper_bound <- m + 3 * sd_x
  x[x < lower_bound | x > upper_bound] <- NA
  return(x)
}

metab4 <- metab4 %>% mutate(across(everything(), replace_outliers))
```

### 8. Filter Metabolites
Exclude metabolites with more than 50% missing values:
```r
metab4 <- metab4[, which(colMeans(!is.na(metab4)) > 0.5)]
```

### 9. Log Transformation and PCA
Perform log transformation and reanalyze PCA:
```r
metab_log <- log(metab2 + 1)
pca_log <- prcomp(metab_log, scale. = TRUE)

ggplot(as.data.frame(pca_log$x), aes(x = PC1, y = PC2)) +
  geom_point() +
  labs(title = "PCA Plot (Log Transformed)", x = "PC1", y = "PC2") +
  theme_dark()
```

### 10. Regression Analysis
Run regression analysis to identify associations between phenotypes and metabolites:
```r
results_matrix <- matrix(0, nrow = ncol(metab4), ncol = 6)
rownames(results_matrix) <- paste("Metabolite", 1:ncol(metab4))
colnames(results_matrix) <- c("Beta_NoCovariate", "SE_NoCovariate", "Pval_NoCovariate",
                              "Beta_WithCovariate", "SE_WithCovariate", "Pval_WithCovariate")

for (i in 1:ncol(metab4)) {
  rs1 <- lm(metab4[, i] ~ phenotype[,"Phenotype1"])
  rs2 <- lm(metab4[, i] ~ phenotype[,"Phenotype1"] + phenotype[,"Phenotype1.1"])
  
  results_matrix[i, ] <- c(
    summary(rs1)$coefficients[2, 1], summary(rs1)$coefficients[2, 2], summary(rs1)$coefficients[2, 4],
    summary(rs2)$coefficients[2, 1], summary(rs2)$coefficients[2, 2], summary(rs2)$coefficients[2, 4]
  )
}
```

### 11. Apply Bonferroni Correction
Identify significant metabolites:
```r
alpha <- 0.05
num_tests <- nrow(results_matrix)
bonferroni_threshold <- alpha / num_tests

significant_metabolites <- results_df[results_df$Pval_NoCovariate < bonferroni_threshold, ]
significant_metabolites <- significant_metabolites[order(significant_metabolites$Pval_NoCovariate), ]
```

### 12. Save Results
Save the significant metabolites to a CSV file:
```r
write.csv(significant_metabolites, "significant_metabolites.csv", row.names = TRUE)
```

---

## Outputs
- PCA visualizations
- Filtered dataset (`metab4`)
- Significant metabolites in `significant_metabolites.csv`

---
```
