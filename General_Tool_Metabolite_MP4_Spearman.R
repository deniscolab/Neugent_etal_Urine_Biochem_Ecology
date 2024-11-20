# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(stats)
library(psych)

# Create output directory
output_dir <- "Metab_Metag_Species"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Read the input CSV files
data_1_T <- read.csv("data_1_T.csv", row.names = 1)
data_2 <- read.csv("Species_relative_abundance.csv", row.names = 1)

# Ensure the samples (rows) are matched between the two datasets
common_samples <- intersect(rownames(data_1_T), rownames(data_2))
data_1_T <- data_1_T[common_samples, ]
data_2 <- data_2[common_samples, ]

# Remove the "class" column from data_1_T if it exists
if ("class" %in% colnames(data_1_T)) {
  data_1_T <- data_1_T %>% select(-class)
}

# Initialize a data frame to store the results
results <- data.frame(Variable1 = character(),
                      Variable2 = character(),
                      Spearman_Correlation = numeric(),
                      p_value = numeric(),
                      adjusted_p_value = numeric(),
                      stringsAsFactors = FALSE)

# Total number of comparisons
total_comparisons <- ncol(data_1_T) * ncol(data_2)
comparison_count <- 0

# Perform Spearman correlation for each pair of columns
for (col1 in colnames(data_1_T)) {
  for (col2 in colnames(data_2)) {
    correlation_test <- cor.test(data_1_T[[col1]], data_2[[col2]], method = "spearman")
    results <- rbind(results, data.frame(Variable1 = col1,
                                         Variable2 = col2,
                                         Spearman_Correlation = correlation_test$estimate,
                                         p_value = correlation_test$p.value))
    # Update progress
    comparison_count <- comparison_count + 1
    progress <- (comparison_count / total_comparisons) * 100
    cat(sprintf("Progress: %.2f%%\n", progress))
  }
}

# Adjust p-values for multiple testing using the Benjamini-Hochberg method
results$adjusted_p_value <- p.adjust(results$p_value, method = "BH")

# Write the results to a CSV file
output_file <- file.path(output_dir, "pairwise_spearman_correlation.csv")
write.csv(results, file = output_file, row.names = FALSE)

# Create a volcano plot
volcano_plot <- ggplot(results, aes(x = Spearman_Correlation, y = -log10(p_value))) +
  geom_point() +
  theme_minimal() +
  labs(title = "Volcano Plot of Spearman Correlation",
       x = "Spearman Correlation",
       y = "-log10(p value)")

# Save the volcano plot
volcano_plot_file <- file.path(output_dir, "volcano_plot.png")
ggsave(volcano_plot_file, plot = volcano_plot)

# Print a message indicating the process is complete
cat("Spearman correlation analysis complete. Results written to", output_dir, "\n")
