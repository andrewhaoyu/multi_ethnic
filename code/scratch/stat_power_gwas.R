# Load required package
library(ggplot2)

# Parameters
cases <- 25000
controls <- 32000
alpha <- 5e-8
maf_list <- c(0.01, 0.05, 0.1)
or_values <- seq(1, 2, by = 0.01)

# Functions to compute NCP and power
compute_ncp <- function(or, maf, n_cases, n_controls) {
  n_eff <- (n_cases * n_controls) / (n_cases + n_controls)
  var_g <- 2 * maf * (1 - maf)
  log_or <- log(or)
  ncp <- n_eff * var_g * (log_or)^2
  return(ncp)
}

compute_power <- function(ncp, alpha) {
  crit <- qchisq(1 - alpha, df = 1)
  power <- 1 - pchisq(crit, df = 1, ncp = ncp)
  return(power)
}

# Generate data frame
power_data <- data.frame()
for (maf in maf_list) {
  for (or in or_values) {
    ncp <- compute_ncp(or, maf, cases, controls)
    power <- compute_power(ncp, alpha)
    power_data <- rbind(power_data, data.frame(OR = or, Power = power, MAF = factor(maf)))
  }
}

# Plot using ggplot2
ggplot(power_data, aes(x = OR, y = Power, color = MAF)) +
  geom_line(size = 1.2) +
  labs(
    title = "Power Curve for Binary Trait GWAS",
    x = "Odds Ratio",
    y = "Statistical Power",
    color = "MAF"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")





# set.seed(42)
# 
# # Simulate genotype given case/control status
# simulate_power_case_control <- function(n_cases, n_controls, maf, or, alpha = 5e-8, n_sim = 1000) {
#   # Genotype frequencies under Hardy-Weinberg Equilibrium
#   pAA <- (1 - maf)^2
#   pAa <- 2 * maf * (1 - maf)
#   paa <- maf^2
#   geno_probs <- c(pAA, pAa, paa)
#   genotypes <- c(0, 1, 2)
# 
#   # Compute disease risk under logistic model
#   log_or <- log(or)
#   logits <- log_or * genotypes
#   probs <- exp(logits) / (1 + exp(logits))  # P(Y=1|G)
# 
#   # Use Bayes’ rule to get P(G | Y=1) and P(G | Y=0)
#   # P(G) = genotype frequencies; P(Y=1|G) = logistic
#   py1 <- sum(probs * geno_probs)
#   py0 <- 1 - py1
#   pG_given_Y1 <- (probs * geno_probs) / py1
#   pG_given_Y0 <- ((1 - probs) * geno_probs) / py0
# 
#   power_count <- 0
# 
#   for (i in 1:n_sim) {
#     # Sample genotypes for cases and controls
#     case_genotypes <- sample(genotypes, size = n_cases, replace = TRUE, prob = pG_given_Y1)
#     control_genotypes <- sample(genotypes, size = n_controls, replace = TRUE, prob = pG_given_Y0)
# 
#     # Combine and run logistic regression
#     genotype_all <- c(case_genotypes, control_genotypes)
#     status_all <- c(rep(1, n_cases), rep(0, n_controls))
# 
#     model <- glm(status_all ~ genotype_all, family = binomial())
#     pval <- summary(model)$coefficients[2, 4]
# 
#     if (!is.na(pval) && pval < alpha) power_count <- power_count + 1
#   }
# 
#   return(power_count / n_sim)
# }
# 
# simulate_power_case_control(
#   n_cases = 25000,
#   n_controls = 32000,
#   maf = 0.01,
#   or = 1.25,
#   alpha = 5e-8,
#   n_sim = 1000
# )
