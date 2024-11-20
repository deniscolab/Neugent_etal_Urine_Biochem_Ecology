# Load necessary library
library(vegan)

# Load the dataset
data <- read.csv("MP4_Relabundance.csv", row.names = 1)

# Function to calculate Hill numbers
calculate_hill_numbers <- function(x) {
  p <- x / sum(x)
  D0 <- sum(x > 0) # 0D - Species Richness
  D1 <- exp(-sum(p * log(p))) # 1D - Shannon Hill Number
  D2 <- 1 / sum(p^2) # 2D - Inverse Simpson Hill Number
  return(c(D0, D1, D2))
}

# Apply the function to each row (sample) and combine the results
hill_numbers <- t(apply(data, 1, calculate_hill_numbers))
colnames(hill_numbers) <- c("0D", "1D", "2D")

# Convert to a data frame
hill_numbers_df <- data.frame(Sample = rownames(data), hill_numbers)

# Create a new directory to save the results
output_dir <- "Hill_Numbers_MP4_Output"
dir.create(output_dir)

# Save the results to a CSV file in the new directory
output_file <- file.path(output_dir, "Hill_Numbers_MP4.csv")
write.csv(hill_numbers_df, output_file, row.names = FALSE)

# Display the results
print(hill_numbers_df)
