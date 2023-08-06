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
PRS = rnorm(N)
data = data.frame(Y,PRS,X_covar)
model1 = lm(Y~PC1+PC2+PC3+PC4+PC5+PC6+
              PC7+PC8+PC9+PC10+sex,data = data)
residual = model1$residuals
model2 = lm(residual~PRS)
summary(model2)$r.square


#for binary traits
Y = rbinom(N,1,0.3)
library(RISCA)
data = data.frame(Y,PRS,X_covar)
Y = rep(c(rep(1,N/2),rep(0,N/2)))
p = 10
X = matrix(rnorm(N*p),N,p)
PRS = rnorm(N)
colnames(X) = paste0("PC",c(1:p))
sex = factor(c(rep("male",N/2),rep("female",N/2)))
X_covar = data.frame(X,sex)

data <- data.frame(Y,PRS,X_covar)



roc_obj = roc.binary(status = "Y",
                     variable = "PRS",
                     confounders = ~PC1+PC2+PC3+PC4+PC5+PC6+
                       PC7+PC8+PC9+PC10+sex,
                     data = data,
                     precision=seq(0.05,0.95, by=0.05))
roc_obj$auc





#code for weighted PRS (Luna 2017, Genetic Epi)
#use PRS matrix from bestEUR PRS
#for linear traits
#Y is outcome
#X are covariates
#two datasets available
#1 represent tuning; 2 represent validation
#PRS is a matrix with four columns
#we only use first column (best EUR) and fourth column (best target)
N = 10000
Y1 = rnorm(N)
p = 10
X1 = matrix(rnorm(N*p),N,p)
colnames(X1) = paste0("PC",c(1:p))
sex1 = factor(c(rep("male",N/2),rep("female",N/2)))
X_covar1 = data.frame(X1,sex = sex1)
PRS_all1 = matrix(rnorm(N*4),N,4)
PRS1 = PRS_all1[,c(1,4),drop=F]
colnames(PRS1) = c("EUR","TAR")
data1 = data.frame(Y = Y1,PRS1,X_covar1)
#data1 is tuning dataset with all 
model1 = lm(Y~PC1+PC2+PC3+PC4+PC5+PC6+
              PC7+PC8+PC9+PC10+sex,data = data1)
residual = model1$residuals
model2 = lm(residual~PRS1)
#get the best weight from tuning dataset
weight = coef(model2)[2:3]
#evaluate the weighted PRS on validation dataset
Y2 = rnorm(N)
X2 = matrix(rnorm(N*p),N,p)
colnames(X2) = paste0("PC",c(1:p))
sex2 = factor(c(rep("male",N/2),rep("female",N/2)))
X_covar2 = data.frame(X2,sex = sex2)
PRS_all2 = matrix(rnorm(N*4),N,4)
PRS2 = PRS_all2[,c(1,4),drop=F]
#Calculate the weighted PRS using weights from tuning data
PRS = PRS2%*%weight
data2 = data.frame(Y = Y2,X_covar2)
model1 = lm(Y~PC1+PC2+PC3+PC4+PC5+PC6+
              PC7+PC8+PC9+PC10+sex,data = data2)
residual = model1$residuals
model2 = lm(residual~PRS)
summary(model2)$r.square


#for binary traits
Y1 = rep(c(rep(1,N/2),rep(0,N/2)))
data1 <- data.frame(Y = Y1,PRS1,X_covar1)

model1 = glm(Y1~PRS1+PC1+PC2+PC3+PC4+PC5+PC6+
               PC7+PC8+PC9+PC10+sex,family= "binomial",
             data=data1)
#ignore the unconvergence message, 
#it's not converged since the simulated data has no signal
weight = coef(model1)[c(2:3)]
#load validation dataset
Y2= rep(c(rep(1,N/2),rep(0,N/2)))
PRS = PRS2%*%weight
data2 = data.frame(Y = Y2,PRS = PRS,X_covar2)
roc_obj = roc.binary(status = "Y",
                     variable = "PRS",
                     confounders = ~PC1+PC2+PC3+PC4+PC5+PC6+
                       PC7+PC8+PC9+PC10+sex,
                     data = data2,
                     precision=seq(0.05,0.95, by=0.05))
roc_obj$auc

#code for PRS-CSx PRS 
#use linear combination PRS matrix with target and eur population
#for linear traits
#Y is outcome
#X are covariates
#two datasets available
#1 represent tuning; 2 represent validation
#PRS is a matrix with eight columns
#1-4 represent PRS based on target population, 5-6 represent PRS based on EUR population
#1 match 5, 2 match 6, 3 match 7, 4 match 8. we have four tuning parameters
N = 10000
Y1 = rnorm(N)
p = 10
X1 = matrix(rnorm(N*p),N,p)
colnames(X1) = paste0("PC",c(1:p))
sex1 = factor(c(rep("male",N/2),rep("female",N/2)))
X_covar1 = data.frame(X1,sex = sex1)
PRS_all1 = matrix(rnorm(N*8),N,8)
data1 = data.frame(Y = Y1,X_covar1)
#record the R2 and regression coefficients
coefficient_mat = matrix(0,4,2)
r2_vec = rep(0,4)
for(k in 1:4){
  PRS1 = PRS_all1[,c(k,k+4),drop=F]
  
  #data1 is tuning dataset with all 
  model1 = lm(Y~PC1+PC2+PC3+PC4+PC5+PC6+
                PC7+PC8+PC9+PC10+sex,data = data1)
  residual = model1$residuals
  model2 = lm(residual~PRS1[,1]+PRS1[,2])
  r2_vec[k] = summary(model2)$r.squared
  coefficient_mat[k,] = coef(model2)[2:3]
}
idx <- which.max(r2_vec)
weight = coefficient_mat[idx,]
#evaluate the weighted PRS on validation dataset
Y2 = rnorm(N)
X2 = matrix(rnorm(N*p),N,p)
colnames(X2) = paste0("PC",c(1:p))
sex2 = factor(c(rep("male",N/2),rep("female",N/2)))
X_covar2 = data.frame(X2,sex = sex2)
PRS_all2 = matrix(rnorm(N*8),N,8)
PRS2 = PRS_all2[,c(idx,idx+4),drop=F]
#Calculate the weighted PRS using weights from tuning data
PRS = PRS2%*%weight
data2 = data.frame(Y = Y2,X_covar2)
model1 = lm(Y~PC1+PC2+PC3+PC4+PC5+PC6+
              PC7+PC8+PC9+PC10+sex,data = data2)
residual = model1$residuals
model2 = lm(residual~PRS)
summary(model2)$r.square




#for binary traits
Y1 = rep(c(rep(1,N/2),rep(0,N/2)))
data1 <- data.frame(Y = Y1,PRS1,X_covar1)
coefficient_mat = matrix(0,4,2)
deviance_vec = rep(0,4)
for(k in 1:4){
  PRS1 = PRS_all1[,c(k,k+4),drop=F]
  
  #data1 is tuning dataset with all 
  model1 = glm(Y1~PRS1+PC1+PC2+PC3+PC4+PC5+PC6+
                 PC7+PC8+PC9+PC10+sex,family= "binomial",
               data=data1)
  deviance_vec[k] = summary(model1)$deviance
  coefficient_mat[k,] = coef(model1)[2:3]
}
idx = which.min(deviance_vec)
weight = coefficient_mat[k,]
PRS2 = PRS_all2[,c(idx,idx+4),drop=F]
#load validation dataset
Y2= rep(c(rep(1,N/2),rep(0,N/2)))
PRS = PRS2%*%weight
data2 = data.frame(Y = Y2,PRS = PRS,X_covar2)
roc_obj = roc.binary(status = "Y",
                     variable = "PRS",
                     confounders = ~PC1+PC2+PC3+PC4+PC5+PC6+
                       PC7+PC8+PC9+PC10+sex,
                     data = data2,
                     precision=seq(0.05,0.95, by=0.05))
roc_obj$auc


#Super learning algorithm
#need both tuning dataset and validation dataset
#for linear traits
#PRS is a matrix from TDLD-EB PRS (972 columns)
library(mvtnorm)
library(dplyr)
library(data.table)
library(caret)
library(SuperLearner)
library(ranger)
library(glmnet)

Y1 = rnorm(N)
p = 10
X1 = matrix(rnorm(N*p),N,p)
colnames(X1) = paste0("PC",c(1:p))
sex1 = factor(c(rep("male",N/2),rep("female",N/2)))
X_covar1 = data.frame(X1,sex = sex1)
#since PRS matrix have correlation
#Here I will use AR1 correlation to simulate the matrix 
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
n.prs = 972
R =ar1_cor(n.prs,0.99)

PRS_all1 = data.frame(rmvnorm(N,mean = rep(0,n.prs),R))
#clean prs.mat
#drop the columns with perfect correlation
mtx = cor(PRS_all1)
drop = findCorrelation(mtx,cutoff = 0.98)
drop = names(PRS_all1)[drop]
prs.mat1 = PRS_all1 %>% 
  select(-all_of(drop))
data1 = data.frame(Y = Y1,X_covar1)
#data1 is tuning dataset with all 
model1 = lm(Y~PC1+PC2+PC3+PC4+PC5+PC6+
              PC7+PC8+PC9+PC10+sex,data = data1)
residual = model1$residuals
SL.libray <- c(
  "SL.glmnet",
  "SL.ridge",
  "SL.nnet"
)
#train super learning model
sl = SuperLearner(Y = residual, X = prs.mat1, family = gaussian(),
                  # For a real analysis we would use V = 10.
                  # V = 3,
                  SL.library = SL.libray)
#load validation dataset
Y2 = rnorm(N)
X2 = matrix(rnorm(N*p),N,p)
colnames(X2) = paste0("PC",c(1:p))
sex2 = factor(c(rep("male",N/2),rep("female",N/2)))
X_covar2 = data.frame(X2,sex = sex2)
PRS_all2 = data.frame(rmvnorm(N,mean = rep(0,n.prs),R))
prs.mat2 = PRS_all2 %>% 
  select(-all_of(drop))
PRS = predict(sl,prs.mat2,onlySL = TRUE)[[1]]
data2 = data.frame(Y = Y2,X_covar2)
model1 = lm(Y~PC1+PC2+PC3+PC4+PC5+PC6+
              PC7+PC8+PC9+PC10+sex,data = data2)
residual = model1$residuals
model2 = lm(residual~PRS)
summary(model2)$r.square




SL.libray <- c(
  "SL.glmnet",
  # "SL.ridge",
  "SL.nnet"
)

#for binary traits
Y1 = rep(c(rep(1,N/2),rep(0,N/2)))
data1 <- data.frame(Y = Y1,PRS1,X_covar1)

sl = SuperLearner(Y = Y1, X = prs.mat1, family = binomial(),
                  # For a real analysis we would use V = 10.
                  # V = 3,
                  SL.library = SL.libray)
#load validation dataset
Y2 = rep(c(rep(1,N/2),rep(0,N/2)))
PRS = predict(sl,prs.mat2,onlySL = TRUE)[[1]]
data2 = data.frame(Y = Y2,PRS = PRS,X_covar2)
roc_obj = roc.binary(status = "Y",
                     variable = "PRS",
                     confounders = ~PC1+PC2+PC3+PC4+PC5+PC6+
                       PC7+PC8+PC9+PC10+sex,
                     data = data2,
                     precision=seq(0.05,0.95, by=0.05))
roc_obj$auc
