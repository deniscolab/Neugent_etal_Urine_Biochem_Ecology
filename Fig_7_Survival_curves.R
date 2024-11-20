# Load necessary libraries
library(survival)
library(survminer)
library(dplyr)

# Create the output directory
if (!dir.exists("time_to_results")) {
  dir.create("time_to_results")
}

# Read the data
data <- read.csv("time_to_data.csv", row.names = 1)

# Extract time to recurrence and recurrence status
time_to_recurrence <- data[, "time.to.recurrence"]
recurrence_status <- data[, "recurrence.status"]

# List of metabolites of interest
metabolites_of_interest <- c("DCA_Bile.Acids", "CE.18.2._Cholesterol.Esters", "CE.16.0._Cholesterol.Esters", "DG.18.2_20.4._Diacylglycerols", "GCDCA_Bile.Acids", "ProBetaine_Aminoacids.Related","PC.ae.C30.2_Glycerophospholipids", "DG.18.2_20.0._Diacylglycerols", "CE.20.5._Cholesterol.Esters", "Cer.d16.1.24.0._Ceramides", "TG.18.3_30.0._Triacylglycerols", "HArg_Aminoacids.Related", "Taurine_Aminoacids.Related", "DG.18.1_20.1._Diacylglycerols", "HexCer.d18.1.18.0._Glycosylceramides", "C3.1_Acylcarnitines", "TG.20.1_30.1._Triacylglycerols", "TG.18.3_38.6._Triacylglycerols", "HexCer.d18.1.23.0._Glycosylceramides", "DG.14.0_18.2._Diacylglycerols", "CE.17.0._Cholesterol.Esters", "HexCer.d18.2.16.0._Glycosylceramides", "DG.14.0_20.0._Diacylglycerols", "TG.22.1_32.5._Triacylglycerols", "PheAlaBetaine_Aminoacids.Related", "lysoPC.a.C24.0_Glycerophospholipids", "lysoPC.a.C26.1_Glycerophospholipids", "DG.22.1_22.2._Diacylglycerols", "PC.aa.C40.1_Glycerophospholipids" )  # Replace with actual metabolite names

# Open a PDF device to save plots
pdf(file = "time_to_results/KM_plots_selected_metabolites.pdf", width = 10, height = 8)

# Loop through each metabolite of interest, dichotomize by median, and create KM plots
for (metabolite in metabolites_of_interest) {
  metabolite_data <- data[[metabolite]]
  
  # Dichotomize by median
  median_value <- median(metabolite_data, na.rm = TRUE)
  dichotomized_data <- ifelse(metabolite_data > median_value, "High", "Low")
  
  # Check for non-missing observations
  if (all(is.na(dichotomized_data))) {
    warning(paste("All values are missing for metabolite:", metabolite))
    next
  }
  
  # Create a survival object
  surv_object <- Surv(time_to_recurrence, recurrence_status)
  
  # Fit the Kaplan-Meier estimator
  fit <- survfit(surv_object ~ dichotomized_data)
  
  # Plot the Kaplan-Meier survival curves
  ggsurv <- ggsurvplot(fit,
                       data = data,
                       ggtheme = theme_minimal(),
                       title = paste("Kaplan-Meier Survival Curve for", metabolite),
                       risk.table = TRUE,
                       conf.int = TRUE,
                       xlab = "Time to Recurrence (days)",
                       ylab = "Survival Probability")
  
  # Print the plot to the PDF
  print(ggsurv)
}

# Close the PDF device
dev.off()

print("Kaplan-Meier plots for selected metabolites saved to time_to_results/KM_plots_selected_metabolites.pdf")
