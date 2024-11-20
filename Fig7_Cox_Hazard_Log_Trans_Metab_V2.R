# Load necessary libraries
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

# Create the output directory
if (!dir.exists("time_to_results")) {
  dir.create("time_to_results")
}

# Create a subdirectory for model summaries
if (!dir.exists("time_to_results/model_summaries")) {
  dir.create("time_to_results/model_summaries")
}

# Read the data
data <- read.csv("time_to_data.csv", row.names = 1)

# Extract time to recurrence and recurrence status
time_to_recurrence <- data[, "time.to.recurrence"]
recurrence_status <- data[, "recurrence.status"]

# Log-transform and standardize the metabolites
metabolite_data <- data[, -c(1, 2)]  # Exclude the first two columns (time to recurrence and recurrence status)
log_transformed_data <- metabolite_data %>%
  mutate(across(everything(), ~ scale(log1p(.))))

# Combine log-transformed data with time to recurrence and recurrence status
log_transformed_data <- cbind(data[, c("time.to.recurrence", "recurrence.status")], log_transformed_data)

# Prepare a data frame to store summary results
results <- data.frame(
  Metabolite = character(),
  Coef = numeric(),
  SE = numeric(),
  Z = numeric(),
  lower_CI = numeric(),
  upper_CI = numeric(),
  PValue = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each metabolite and fit the Cox model
for (metabolite in colnames(log_transformed_data)[-c(1, 2)]) {
  cox_model <- coxph(Surv(time_to_recurrence, recurrence_status) ~ log_transformed_data[[metabolite]], data = log_transformed_data)
  summary_model <- summary(cox_model)
  
  # Extract the coefficient, standard error, z-value, confidence intervals, and p-value
  coef <- summary_model$coef[1, "coef"]
  se_coef <- summary_model$coef[1, "se(coef)"]
  z_value <- summary_model$coef[1, "z"]
  ci_lower <- summary_model$conf.int[1, "lower .95"]
  ci_upper <- summary_model$conf.int[1, "upper .95"]
  p_value <- summary_model$coef[1, "Pr(>|z|)"]
  
  # Store the results in the results data frame
  results <- rbind(results, data.frame(
    Metabolite = metabolite,
    Coef = coef,
    SE = se_coef,
    Z = z_value,
    lower_CI = ci_lower,
    upper_CI = ci_upper,
    PValue = p_value,
    stringsAsFactors = FALSE
  ))
  
  # Save the full summary model as a CSV
  summary_df <- as.data.frame(summary_model$coefficients)
  summary_filename <- paste0("time_to_results/model_summaries/summary_", metabolite, ".csv")
  write.csv(summary_df, summary_filename, row.names = TRUE)
}

# Write the summary results to a CSV file
write.csv(results, "time_to_results/cox_model_results_log_transformed_corrected.csv", row.names = FALSE)

# Create the volcano plot
results$NegLogPValue = -log10(results$PValue)
volcano_plot <- ggplot(results, aes(x = Coef, y = NegLogPValue)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  labs(title = "Volcano Plot of Metabolites",
       x = "Log Hazard Ratio (Coefficient)",
       y = "-Log10(P-value)") +
  theme_minimal()

# Save the volcano plot
ggsave(filename = "time_to_results/Volcano_Plot_log_transformed_corrected.pdf", plot = volcano_plot, width = 10, height = 8)

print("Cox model analysis complete. Results saved to time_to_results/cox_model_results_log_transformed_corrected.csv and Volcano_Plot_log_transformed_corrected.pdf")
print("Full summaries saved in time_to_results/model_summaries/")
