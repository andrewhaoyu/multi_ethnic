#for linear traits
#Y is outcome
#X are covariates
#PRS is one column
N = 10000
Y = rnorm(N)
p = 10
X = matrix(rnorm(N*p),N,p)
colnames(X) = paste0("PC",c(1:p))
sex = factor(c(rep("male",N/2),rep("female",N/2)))
X_covar = data.frame(X,sex)
data = data.frame(Y,PRS,X_covar)
PRS = rnorm(N)
model1 = lm(Y~PC1+PC2+PC3+PC4+PC5+PC6+
              PC7+PC8+PC9+PC10+sex,data = data)
residual = model1$residuals
model2 = lm(residual~PRS)
summary(model2)$r.square


#for binary traits
Y = rbinom(N,1,0.3)
library(RISCA)
data = data.frame(Y,PRS,X_covar)
Y = rep(c(rep("yes",N/2),rep("no",N/2)))
data <- data.frame(Y,PRS,X_covar)



roc_obj = roc.binary(status = "Y",
                     variable = "PRS",
                     confounders = PC1+PC2+PC3+PC4+PC5+PC6+
                       PC7+PC8+PC9+PC10+sex,
                     data = data,
                     precision=seq(0.05,0.95, by=0.05))
roc_obj$auc
