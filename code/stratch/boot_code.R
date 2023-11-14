# Load necessary library
library(lmtest)

# Generate synthetic data
set.seed(123)
n <- 100 # Number of observations
X <- rnorm(n, 50, 10)
Y <- 0.1 * X + rnorm(n, 0, 5) # Y is linearly related to X
covariate1 <- rnorm(n, 0, 1)
covariate2 <- rnorm(n, 0, 1)
data <- data.frame(X, Y, covariate1, covariate2)

# Fit the initial linear model
model <- lm(Y ~ X + covariate1 + covariate2, data = data)
initial_slope <- coef(model)["X"]

# Permutation procedure
n_permutations <- 10000
permuted_slopes <- numeric(n_permutations)

for(i in 1:n_permutations) {
  permuted_Y <- sample(data$Y) # Shuffle Y values
  permuted_data <- data
  permuted_data$Y <- permuted_Y
  
  permuted_model <- lm(Y ~ X + covariate1 + covariate2, data = permuted_data)
  permuted_slopes[i] <- coef(permuted_model)["X"]
}

# Calculate the two-sided p-value
p_value <- mean(abs(permuted_slopes) >= abs(initial_slope))

# Output the results
cat("Permutation P-Value:", p_value, "\n")
cat("Model-based P-Value:", coefficients(summary(model))[2,4], "\n")
