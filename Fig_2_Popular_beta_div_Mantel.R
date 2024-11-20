# Load necessary libraries
library(vegan)  # For distance calculations and Mantel test
library(dplyr)  # For data manipulation

# Load the datasets
metabolites <- read.csv("BC_data_full_cohort_working.csv", row.names = 1)
taxa <- read.csv("Taxa_MP4_Relabundance.csv", row.names = 1)

# Check if sample names match in both datasets
if (!all(rownames(metabolites) == rownames(taxa))) {
  stop("Sample names do not match between datasets.")
}

# List of popular beta diversity metrics, excluding Minkowski
beta_metrics <- c("bray", "jaccard", "euclidean", "canberra", "kulczynski", "gower", "manhattan")

# Initialize a data frame to store the results
results <- data.frame(metric_metabolites = character(), metric_taxa = character(), 
                      mantel_statistic = numeric(), p_value = numeric())

# Iterate over all combinations of beta diversity metrics
for (metric_metabolites in beta_metrics) {
  for (metric_taxa in beta_metrics) {
    
    # Calculate beta diversity distance for metabolites dataset
    metabolite_distance <- vegdist(metabolites, method = metric_metabolites)
    
    # Calculate beta diversity distance for taxa dataset
    taxa_distance <- vegdist(taxa, method = metric_taxa)
    
    # Perform the Spearman Mantel test between the two distance matrices
    mantel_test_result <- mantel(as.dist(metabolite_distance), as.dist(taxa_distance), method = "spearman")
    
    # Store the metric names, Mantel test statistic, and p-value in the results data frame
    results <- rbind(results, data.frame(
      metric_metabolites = metric_metabolites,
      metric_taxa = metric_taxa,
      mantel_statistic = mantel_test_result$statistic,
      p_value = mantel_test_result$signif
    ))
  }
}

# Output the results to a CSV file
write.csv(results, "mantel_test_results.csv", row.names = FALSE)

# Print a message indicating the completion of the process
print("Mantel test results saved to 'mantel_test_results.csv'")
