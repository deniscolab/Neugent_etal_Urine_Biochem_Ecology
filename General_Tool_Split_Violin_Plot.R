library(vioplot)

# Load the data
data <- read.csv("Alpha_Diversity_Full_Cohort.csv")

# Define custom colors for each group (lighter for left, darker for right)
colors_left <- c("lightblue", "#DDA0DD", "lightcoral")  # Light shades for left
colors_right <- c("deepskyblue", "purple", "red")       # Darker shades for right

# Split the data by class
classes <- unique(data$class)

# Save as PDF for Adobe Illustrator
pdf("split_violin_plot.pdf")

# Set up the plot area
plot(1, type = "n", xlim = c(0.5, length(classes) + 0.5), 
     ylim = range(c(data$Biochemical, data$Taxonomic)),
     xaxt = "n", xlab = "Class", ylab = "Alpha Diversity", main = "Split Violin Plot")

# Custom x-axis labels
axis(1, at = 1:length(classes), labels = classes)

# Loop through each class and create split violins
for (i in 1:length(classes)) {
  bio_data <- data[data$class == classes[i], "Biochemical"]
  tax_data <- data[data$class == classes[i], "Taxonomic"]
  
  # Plot left and right split violins with custom colors
  vioplot(bio_data, at = i, side = "left", col = colors_left[i], add = TRUE, plotCentre = "line")
  vioplot(tax_data, at = i, side = "right", col = colors_right[i], add = TRUE, plotCentre = "line")
}

# Close the PDF device to save the file
dev.off()
