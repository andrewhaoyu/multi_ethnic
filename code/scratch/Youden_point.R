# Install and load the pROC package
#install.packages("pROC")
library(pROC)

# Generate some example data
set.seed(123)
response <- sample(c(0, 1), 100, replace = TRUE)
predictor <- rnorm(100)

# Create a ROC curve
roc_curve <- roc(response, predictor)

# Find the threshold that maximizes the Youden's index
best_threshold <- coords(roc_curve, "best")

print(best_threshold)

#Interpretation 1: sensitivity
best_threshold$sensitivity


#Interpretation 2: classified true cases among all cases
predicted_classes <- ifelse(predictor >= best_threshold$threshold, 1, 0)
accuracy <- sum(predicted_classes == response) / length(response)
print(accuracy)
