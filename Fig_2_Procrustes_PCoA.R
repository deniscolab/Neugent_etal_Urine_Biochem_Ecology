# Load necessary libraries
library(vegan)  # For distance calculations and Procrustes analysis
library(ggplot2)  # For plotting

# Load the datasets
metabolites <- read.csv("BC_data_full_cohort_working.csv", row.names = 1)
taxa <- read.csv("Taxa_MP4_Relabundance.csv", row.names = 1)

# Check if sample names match in both datasets
if (!all(rownames(metabolites) == rownames(taxa))) {
  stop("Sample names do not match between datasets.")
}

# Create a new directory for output (if it doesn't exist)
output_dir <- "Procrustes_Results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# List of popular beta diversity metrics, excluding Minkowski
beta_metrics <- c("bray", "jaccard", "euclidean", "canberra", "kulczynski", "gower", "manhattan")

# Initialize a data frame to store the summary (Procrustes sum of squares, m_squared, protest p-value)
summary_results <- data.frame(metric_metabolites = character(), 
                              metric_taxa = character(), 
                              m_squared = numeric(),
                              protest_p_value = numeric())

# Iterate over all combinations of beta diversity metrics
for (metric_metabolites in beta_metrics) {
  for (metric_taxa in beta_metrics) {
    
    # Calculate beta diversity distance for metabolites dataset
    metabolite_distance <- vegdist(metabolites, method = metric_metabolites)
    
    # Calculate beta diversity distance for taxa dataset
    taxa_distance <- vegdist(taxa, method = metric_taxa)
    
    # Perform Principal Coordinates Analysis (PCoA) for ordination
    pcoa_metabolites <- cmdscale(metabolite_distance, k = 4)  # k=4 dimensions
    pcoa_taxa <- cmdscale(taxa_distance, k = 4)  # k=4 dimensions
    
    # Perform the Procrustes analysis between the two ordination results
    procrustes_result <- procrustes(pcoa_metabolites, pcoa_taxa)
    
    # Perform Procrustes permutation test (protest) to assess significance
    protest_result <- protest(pcoa_metabolites, pcoa_taxa, permutations = 999)
    
    # Extract the Procrustes sum of squares (m_squared) and p-value from protest
    m_squared <- procrustes_result$ss
    protest_p_value <- protest_result$signif
    
    # Append the m_squared value and protest p-value to the summary data frame
    summary_results <- rbind(summary_results, data.frame(
      metric_metabolites = metric_metabolites,
      metric_taxa = metric_taxa,
      m_squared = m_squared,
      protest_p_value = protest_p_value
    ))
    
    # Save detailed Procrustes results (rotation matrix, residuals) to individual CSV files
    procrustes_output_dir <- file.path(output_dir, paste0(metric_metabolites, "_vs_", metric_taxa))
    if (!dir.exists(procrustes_output_dir)) {
      dir.create(procrustes_output_dir)
    }
    
    # Save rotation matrix to CSV
    write.csv(procrustes_result$rotation, file.path(procrustes_output_dir, "procrustes_rotation.csv"), row.names = FALSE)
    
    # Save residuals to CSV
    write.csv(procrustes_result$residuals, file.path(procrustes_output_dir, "procrustes_residuals.csv"), row.names = FALSE)
    
    # --- Add the Procrustes Alignment Plot as an Editable PDF ---
    # Create the PDF output file for Procrustes alignment plot (set 2x2 inches plot area)
    pdf(file = file.path(procrustes_output_dir, "Procrustes_Alignment.pdf"), width = 2, height = 2)
    
    # Set margins to ensure the graph itself is 2x2 inches
    par(mar = c(1, 1, 1, 1))  # Adjust margins to maximize plotting area
    
    # Plot the default Procrustes alignment
    plot(procrustes_result, main = paste("Procrustes Alignment:", metric_metabolites, "vs", metric_taxa))
    
    # Close the PDF device to finalize the file
    dev.off()
  }
}

# Save the summary results (with protest p-value) to a CSV file in the new directory
write.csv(summary_results, file.path(output_dir, "procrustes_summary_with_p_value.csv"), row.names = FALSE)

# Print a message indicating the completion of the process
print(paste("Procrustes results, alignment plots, and summary with p-values saved in", output_dir))
