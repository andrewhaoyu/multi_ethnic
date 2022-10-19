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
data <- data.frame(Y,PRS,X_covar)



roc_obj = roc.binary(status = "Y",
                     variable = "PRS",
                     confounders = ~PC1+PC2+PC3+PC4+PC5+PC6+
                       PC7+PC8+PC9+PC10+sex,
                     data = data,
                     precision=seq(0.05,0.95, by=0.05))
roc_obj$auc


#24 columns in the three methods prs data
#the 1 column is the PRS for XPASS
#the 2-4 columns is for POLYPRED+
#the 5-24 columns are for PRS-CSx


#two datasets available
#1 represent tuning; 2 represent validation
#PRS is a matrix with 24 columns
PRS_all1 = matrix(rnorm(N*24),N,24)
PRS_all2 = matrix(rnorm(N*24),N,24)

#XPASS method (Cai, AJHG, 2021)
#take the PRS for XPASS and directly evaluates it on validation data
#for continuous traits
Y2 = rnorm(N)
X2 = matrix(rnorm(N*p),N,p)
colnames(X2) = paste0("PC",c(1:p))
sex2 = factor(c(rep("male",N/2),rep("female",N/2)))
X_covar2 = data.frame(X2,sex = sex2)
PRS = PRS_all2[,1,drop=F]
#Calculate the weighted PRS using weights from tuning data
data2 = data.frame(Y = Y2,X_covar2)
model1 = lm(Y~PC1+PC2+PC3+PC4+PC5+PC6+
              PC7+PC8+PC9+PC10+sex, data = data2)
residual = model1$residuals
model2 = lm(residual~PRS)
summary(model2)$r.square

#for binary traits
#load validation dataset
Y2= rep(c(rep(1,N/2),rep(0,N/2)))

data2 = data.frame(Y = Y2,PRS = PRS,X_covar2)
roc_obj = roc.binary(status = "Y",
                     variable = "PRS",
                     confounders = ~PC1+PC2+PC3+PC4+PC5+PC6+
                       PC7+PC8+PC9+PC10+sex,
                     data = data2,
                     precision=seq(0.05,0.95, by=0.05))
roc_obj$auc




#PolyPred+ method (Weissbrod, NG, 2022)
#column 3-5 is for 
#Estimates the weight for the three PRSs and evaluate on the validation
#for continuous traits
Y2 = rnorm(N)
PRS1 = PRS_all1[,c(2:4),drop=F]
#if a PRS doesn't converge, I put all the coefficients as 0 for formating
col_sum = colSums(PRS1)
#find the columns that converges
idx <- which(col_sum!=0)
PRS1 = PRS1[, idx, drop = F]

data1 = data.frame(Y = Y1,PRS1,X_covar1)
#data1 is tuning dataset with all 
model1 = lm(Y~PC1+PC2+PC3+PC4+PC5+PC6+
              PC7+PC8+PC9+PC10+sex,data = data1)
residual = model1$residuals
model2 = lm(residual~PRS1)
#get the best weight from tuning dataset
weight = coef(model2)[2:(ncol(PRS1)+1)]
#evaluate the weighted PRS on validation dataset
PRS2 = as.matrix(PRS_all2[,c(2:4),drop=F])
PRS2 = PRS2[, idx, drop = F]
#Calculate the weighted PRS using weights from tuning data
PRS = PRS2%*%weight
data2 = data.frame(Y = Y2,X_covar2)
model1 = lm(Y~PC1+PC2+PC3+PC4+PC5+PC6+
              PC7+PC8+PC9+PC10+sex,data = data2)
residual = model1$residuals
model2 = lm(residual~PRS)
summary(model2)$r.square


#for binary traits
#load validation dataset
Y2= rep(c(rep(1,N/2),rep(0,N/2)))

data2 = data.frame(Y = Y2,PRS = PRS,X_covar2)
roc_obj = roc.binary(status = "Y",
                     variable = "PRS",
                     confounders = ~PC1+PC2+PC3+PC4+PC5+PC6+
                       PC7+PC8+PC9+PC10+sex,
                     data = data2,
                     precision=seq(0.05,0.95, by=0.05))
roc_obj$auc




#code for PRS-CSx PRS (five ancestries)
#column 5-24 are for PRS-csx
#PRS-CSx has a tuning parameter phi with four values
#we need to find the optimal phi 
#under each phi, we will linear weight the PRS of five ancestries
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
PRS_all1 = matrix(rnorm(N*24),N,24)
data1 = data.frame(Y = Y1,X_covar1)
#record the R2 and regression coefficients
coefficient_mat = matrix(0,4,5)
r2_vec = rep(0,4)
PRS_csx_all1 = PRS_all1[,c(5:24)]
for(k in 1:4){
  PRS1 = PRS_csx_all1[,5*(k-1)+c(1:5),drop=F]
  
  #data1 is tuning dataset with all 
  model1 = lm(Y~PC1+PC2+PC3+PC4+PC5+PC6+
                PC7+PC8+PC9+PC10+sex,data = data1)
  residual = model1$residuals
  model2 = lm(residual~PRS1)
  r2_vec[k] = summary(model2)$r.squared
  coefficient_mat[k,] = coef(model2)[2:(ncol(PRS1)+1)]
}
idx <- which.max(r2_vec)
weight = coefficient_mat[idx,]
#evaluate the weighted PRS on validation dataset
Y2 = rnorm(N)
PRS_csx_all2 = PRS_all2[,c(5:24)]
PRS2 = as.matrix(PRS_csx_all2[,5*(idx-1)+c(1:5)])
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
coefficient_mat = matrix(0,4,5)
deviance_vec = rep(0,4)
PRS_csx_all1 = PRS_all1[,c(5:24)]
for(k in 1:4){
  PRS1 = PRS_csx_all1[,5*(k-1)+c(1:5),drop=F]
  
  #data1 is tuning dataset with all 
  model1 = glm(Y1~PRS1+PC1+PC2+PC3+PC4+PC5+PC6+
                 PC7+PC8+PC9+PC10+sex,family= "binomial",
               data=data1)
  deviance_vec[k] = summary(model1)$deviance
  coefficient_mat[k,] = coef(model1)[2:(ncol(PRS1)+1)]
}
idx = which.min(deviance_vec)
weight = coefficient_mat[k,]
PRS_csx_all2 = PRS_all2[,c(5:24)]
PRS2 = as.matrix(PRS_csx_all2[,5*(idx-1)+c(1:5)])
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


