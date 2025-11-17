# Fixed-effect meta-analysis in R

# Input: vectors of betas and standard errors
beta <- c(1.7, -0.1)
se <- c(0.3, 0.94)

# Step 1: compute weights
weights <- 1 / (se^2)

# Step 2: fixed-effect meta-analysis estimate
beta_meta <- sum(weights * beta) / sum(weights)

# Step 3: standard error of meta-analysis
se_meta <- sqrt(1 / sum(weights))

# Step 4: Z statistic and p-value
z_meta <- beta_meta / se_meta
p_meta <- 2 * pnorm(-abs(z_meta))

# Output results
cat("Fixed-effect meta-analysis results:\n")
cat("Combined Beta:", round(beta_meta, 4), "\n")
cat("Combined SE:", round(se_meta, 4), "\n")
cat("Z-score:", round(z_meta, 4), "\n")
cat("P-value:", format(p_meta, scientific = TRUE), "\n")
