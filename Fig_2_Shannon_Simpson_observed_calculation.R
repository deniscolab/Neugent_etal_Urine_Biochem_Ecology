# Load necessary libraries
library(vegan)
library(dplyr)

# Load the data
data <- read.csv("BC_Data_G1G3.csv", row.names = 1)

# Calculate alpha diversity metrics
shannon <- diversity(data, index = "shannon")
simpson <- diversity(data, index = "simpson")
observed <- rowSums(data > 0)

# Combine the results into a data frame
alpha_diversity <- data.frame(Sample = rownames(data),
                              Shannon = shannon,
                              Simpson = simpson,
                              Observed = observed)

# Write the results to a new CSV file
write.csv(alpha_diversity, "alpha_diversity_metrics.csv", row.names = FALSE)

# Print message indicating completion
cat("Alpha diversity metrics calculated and saved to alpha_diversity_metrics.csv")
