library(boot)
library(RISCA)
# define a function that calculates the AUC
calc_auc <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  roc_obj = roc.binary(status = "Y",
                       variable = "PRS",
                       confounders = ~PC1+PC2+PC3+PC4+PC5+PC6+
                         PC7+PC8+PC9+PC10+sex,
                       data = d,
                       precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}

#data generation process
N = 100 # example number of observations
p = 10
Y = rep(c(rep(1,N/2),rep(0,N/2)))
X = matrix(rnorm(N*p),N,p)
PRS = rnorm(N)
colnames(X) = paste0("PC",c(1:p))
sex = factor(c(rep("male",N/2),rep("female",N/2)))
X_covar = data.frame(X,sex)

data <- data.frame(Y,PRS,X_covar)

# perform bootstrap resampling
set.seed(123)  # for reproducibility
boot.results <- boot(data=data, statistic=calc_auc, R=2000)  # R is the number of bootstrap replicates

# calculate 95% CI with BCA method
boot.ci(boot.out=boot.results, conf=0.95, type="bca")
