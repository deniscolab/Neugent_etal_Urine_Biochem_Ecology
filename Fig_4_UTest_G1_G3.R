# Load necessary libraries
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(dplyr)
library(tidyr)
library(ggplot2)

# Load data
data <- read.csv("data_1_T.csv", header = TRUE, row.names = 1)

# Extract class information
classes <- as.factor(data[, 1])
data <- data[, -1]

# Convert data to numeric
data <- as.data.frame(lapply(data, as.numeric))

# Subset data for classes 1 and 3
class_1_data <- data[which(classes == 1), ]
class_3_data <- data[which(classes == 3), ]

# Function for non-parametric differential analysis
perform_analysis <- function(feature) {
  # Median values for each class
  median_class_1 <- median(class_1_data[[feature]], na.rm = TRUE)
  median_class_3 <- median(class_3_data[[feature]], na.rm = TRUE)
  
  # Median fold change
  median_fc <- median_class_1 / median_class_3
  
  # Wilcoxon rank sum test
  test <- wilcox.test(class_1_data[[feature]], class_3_data[[feature]])
  
  # Z value (effect size)
  z_value <- qnorm(test$p.value / 2) * sign(median_class_1 - median_class_3)
  
  # Mean rank difference
  ranks <- rank(c(class_1_data[[feature]], class_3_data[[feature]]))
  mean_rank_diff <- mean(ranks[1:nrow(class_1_data)]) - mean(ranks[(nrow(class_1_data) + 1):(nrow(class_1_data) + nrow(class_3_data))])
  
  # Determine enrichment
  enrichment <- ifelse(median_class_1 > median_class_3, "Class 1", "Class 3")
  
  # Return results
  return(c(median_fc, z_value, mean_rank_diff, test$p.value, enrichment))
}

# Apply analysis function to each feature
results <- t(sapply(colnames(data), perform_analysis))
colnames(results) <- c("Median_Fold_Change", "Z_Value", "Mean_Rank_Difference", "P_Value", "Enrichment")

# Adjust p-values for multiple comparisons
results <- as.data.frame(results)
results$Adjusted_P_Value <- p.adjust(as.numeric(results$P_Value), method = "fdr")

# Reorder columns to include adjusted P value before enrichment
results <- results %>%
  select(Median_Fold_Change, Z_Value, Mean_Rank_Difference, P_Value, Adjusted_P_Value, Enrichment)

# Output results
write.csv(results, "differential_analysis_results.csv", row.names = TRUE)

# Display results
print(results)
