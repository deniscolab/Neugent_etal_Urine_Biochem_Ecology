# Load necessary libraries
library(tidyverse)
library(psych)
library(ggpubr)
library(patchwork)

# Set file paths
file1_path <- "MP4_genera_abundance_mclr_G1G3.csv"
file2_path <- "Biocrates_Discriminating_lipids_data_T_G1G3.csv"

# Read data
data1 <- read.csv(file1_path, row.names = 1)
data2 <- read.csv(file2_path, row.names = 1)

# Ensure row names match
common_samples <- intersect(rownames(data1), rownames(data2))
data1 <- data1[common_samples, ]
data2 <- data2[common_samples, ]

# Output directory
output_dir <- "Spearman_Correlation_Results_Genera_V3"
dir.create(output_dir, showWarnings = FALSE)

# Initialize result dataframe
results <- data.frame(
  Feature1 = character(),
  Feature2 = character(),
  Spearman_R = numeric(),
  P_value = numeric(),
  Padj_value = numeric(),
  stringsAsFactors = FALSE
)

# Initialize list for significant plots
significant_plots <- list()

# Calculate Spearman correlation
total_combinations <- ncol(data1) * ncol(data2)
combination_counter <- 0

for (i in seq_along(colnames(data1))) {
  for (j in seq_along(colnames(data2))) {
    combination_counter <- combination_counter + 1
    cat("Processing combination", combination_counter, "of", total_combinations, "\n")
    
    feature1 <- colnames(data1)[i]
    feature2 <- colnames(data2)[j]
    
    # Perform Spearman correlation
    corr_test <- corr.test(data1[[feature1]], data2[[feature2]], method = "spearman")
    spearman_r <- corr_test$r
    p_value <- corr_test$p
    
    # Append results if significant and greater than 0.3
    if (p_value < 0.05 && abs(spearman_r) > 0.3) {
      results <- rbind(results, data.frame(
        Feature1 = feature1,
        Feature2 = feature2,
        Spearman_R = spearman_r,
        P_value = p_value,
        stringsAsFactors = FALSE
      ))
      
      # Save significant plot
      plot_data <- data.frame(x = data1[[feature1]], y = data2[[feature2]])
      colnames(plot_data) <- c(feature1, feature2)
      
      plot <- ggscatter(plot_data,
                        x = feature1,
                        y = feature2,
                        add = "reg.line",
                        conf.int = TRUE,
                        cor.coef = TRUE,
                        cor.method = "spearman",
                        xlab = feature1,
                        ylab = feature2,
                        title = paste("Spearman R:", round(spearman_r, 3), "P:", p_value)
      )
      significant_plots[[length(significant_plots) + 1]] <- plot
    }
  }
}

# Adjust p-values
results$Padj_value <- p.adjust(results$P_value, method = "fdr")

# Write results to CSV
write.csv(results, file.path(output_dir, "spearman_correlation_results.csv"), row.names = FALSE)

# Save significant plots to PDF
pdf(file.path(output_dir, "significant_scatter_plots.pdf"))
for (plot in significant_plots) {
  print(plot)
}
dev.off()

cat("Analysis complete. Results saved to", output_dir, "\n")
