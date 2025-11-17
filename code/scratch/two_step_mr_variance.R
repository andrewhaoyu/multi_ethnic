# Input MR estimates and SEs
beta1 <- 0.2     # AS → IL17
se1 <- 0.05

beta2 <- 0.3     # IL17 → Pancreatic Cancer
se2 <- 0.07

# Indirect effect
indirect <- beta1 * beta2

# Standard error via delta method
se_indirect <- sqrt((beta2^2 * se1^2) + (beta1^2 * se2^2))

# P value
z <- indirect/se_indirect
p_value <- 2 * pnorm(-abs(z), lower.tail = T)

# 95% CI
ci_lower <- indirect - 1.96 * se_indirect
ci_upper <- indirect + 1.96 * se_indirect

# Output
cat("Indirect effect:", round(indirect, 4), "\n")
cat("95% CI:", round(ci_lower, 4), "-", round(ci_upper, 4), "\n")
