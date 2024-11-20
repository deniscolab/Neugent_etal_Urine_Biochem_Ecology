# Load necessary libraries if not already loaded
if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")

library(glmnet)
library(caret)
library(pROC)

# Set seed for reproducibility
set.seed(12345)

# Load data
data <- read.csv("summed_lipids.csv")

# Subset data for classes 1 and 3
data_1_3 <- data[data$Class %in% c(1, 3), ]

# Recode Class as 0 and 1
data_1_3$Class <- ifelse(data_1_3$Class == 1, 0, 1)

# Define features (excluding "Sample", "Class", and "Summed_Lipids" columns)
features <- colnames(data_1_3)[!colnames(data_1_3) %in% c("Sample", "Class", "Summed_Lipids")]

# Convert to matrix format for glmnet
x <- as.matrix(data_1_3[, features])
y <- data_1_3$Class

# Set a more detailed grid for alpha values
alpha_values <- seq(0, 1, by = 0.1)

# Perform repeated cross-validated Elastic Net to find the best alpha and lambda
cv_fit <- lapply(alpha_values, function(alpha) {
  cv.glmnet(x, y, family = "binomial", alpha = alpha, nfolds = 10, type.measure = "auc")
})

# Find the best combination of alpha and lambda
best_fit <- cv_fit[[which.max(sapply(cv_fit, function(fit) max(fit$cvm)))]]
best_alpha <- alpha_values[which.max(sapply(cv_fit, function(fit) max(fit$cvm)))]
best_lambda <- best_fit$lambda.min

print(paste("Best alpha:", best_alpha))
print(paste("Best lambda:", best_lambda))

# Extract the coefficients for the best model
best_model <- glmnet(x, y, family = "binomial", alpha = best_alpha, lambda = best_lambda)
coefficients <- coef(best_model)
optimal_features <- rownames(coefficients)[which(coefficients != 0)]
optimal_features <- optimal_features[optimal_features != "(Intercept)"]

print(optimal_features)

# Calculate variable importance (absolute value of coefficients)
importance <- abs(as.vector(coefficients[which(coefficients != 0)]))
importance <- importance[-1]  # Remove intercept

variable_importance <- data.frame(Feature = optimal_features, Importance = importance)

# Save the optimal features and their importance to a new CSV file
output_dir <- "Elastic_Net"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
write.csv(variable_importance, file.path(output_dir, "optimal_features_importance.csv"), row.names = FALSE)

# Initialize vectors to store LOOCV results
loocv_probabilities <- numeric(nrow(data_1_3))
loocv_classes <- numeric(nrow(data_1_3))

# Perform LOOCV using the optimal Elastic Net model
for (i in 1:nrow(data_1_3)) {
  # Split data into training and test sets
  train_indices <- setdiff(1:nrow(data_1_3), i)
  x_train <- x[train_indices, ]
  y_train <- y[train_indices]
  x_test <- x[i, , drop = FALSE]
  
  # Fit Elastic Net model on training data
  elastic_net_model <- glmnet(x_train, y_train, family = "binomial", alpha = best_alpha, lambda = best_lambda)
  
  # Predict probability for the left-out observation
  loocv_probabilities[i] <- predict(elastic_net_model, newx = x_test, type = "response")
}

# Determine predicted classes based on probability threshold of 0.5
loocv_classes <- ifelse(loocv_probabilities > 0.5, 1, 0)
data_1_3$LOOCV_Predicted_Class <- factor(loocv_classes, levels = c(0, 1))

# Confusion Matrix for LOOCV
conf_matrix_loocv <- confusionMatrix(data_1_3$LOOCV_Predicted_Class, factor(data_1_3$Class, levels = c(0, 1)))
print(conf_matrix_loocv)

# ROC Curve and AUC for LOOCV
roc_curve_loocv <- roc(data_1_3$Class, loocv_probabilities, levels = c(0, 1), direction = "<")
auc_value_loocv <- auc(roc_curve_loocv)
print(auc_value_loocv)

# Save ROC Curve coordinates to CSV
roc_coords <- data.frame(Specificity = roc_curve_loocv$specificities, Sensitivity = roc_curve_loocv$sensitivities)
write.csv(roc_coords, file.path(output_dir, "roc_curve_coordinates.csv"), row.names = FALSE)

# Perform logistic regression on the full dataset using the optimal features
full_logistic_model <- glm(Class ~ ., data = data_1_3[, c(optimal_features, "Class")], family = binomial)

# Perform Likelihood Ratio Test
null_model <- glm(Class ~ 1, data = data_1_3[, c(optimal_features, "Class")], family = binomial)
lrt <- anova(null_model, full_logistic_model, test = "LRT")
lrt_p_value <- lrt$`Pr(>Chi)`[2]

# Save the LRT results
lrt_results <- data.frame(Metric = "Likelihood Ratio Test P-Value", Value = lrt_p_value)
write.csv(lrt_results, file.path(output_dir, "lrt_results.csv"), row.names = FALSE)

# Save full logistic regression model coefficients
logistic_model_coefficients <- summary(full_logistic_model)$coefficients
write.csv(logistic_model_coefficients, file.path(output_dir, "logistic_model_coefficients.csv"), row.names = TRUE)

# Save the results to CSV files
# Save ROC Curve data to CSV
roc_curve_data <- data.frame(True_Class = data_1_3$Class, Predicted_Probability = loocv_probabilities)
write.csv(roc_curve_data, file.path(output_dir, "roc_curve_data.csv"), row.names = FALSE)

# Save confusion matrix data to CSV
conf_matrix_data <- as.data.frame(conf_matrix_loocv$table)
write.csv(conf_matrix_data, file.path(output_dir, "confusion_matrix.csv"), row.names = FALSE)

# Save model predictions and actuals
predictions <- data.frame(Actual_Class = data_1_3$Class, Predicted_Class = data_1_3$LOOCV_Predicted_Class)
write.csv(predictions, file.path(output_dir, "predictions.csv"), row.names = FALSE)
