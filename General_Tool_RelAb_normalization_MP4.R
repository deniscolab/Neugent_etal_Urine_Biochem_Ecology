# Load necessary libraries
library(dplyr)

# Read the CSV file
input_file <- "MP4_Genera_Counts.csv"
data <- read.csv(input_file, row.names = 1)

# Normalize each cell to the sum of its row
normalized_data <- data %>%
  mutate_all(~ ./rowSums(data))

# Write the normalized data to a new CSV file
output_file <- "MP4_genera_abundance.csv"
write.csv(normalized_data, file = output_file, row.names = TRUE)

# Print a message indicating the process is complete
cat("Normalization complete. Output written to", output_file, "\n")
