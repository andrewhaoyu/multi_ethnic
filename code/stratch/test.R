#test pcakge
install.packages("SuperLearner")
install.packages(c("caret", "glmnet", "randomForest", "ggplot2", "RhpcBLASctl"))
install.packages("xgboost", repos=c("http://dmlc.ml/drat/", getOption("repos")), type="source")
data(Boston, package = "MASS")
library(ranger)
colSums(is.na(Boston))
outcome = Boston$medv
data = subset(Boston, select = -medv)
str(data)
dim(data)
set.seed(1)
train_obs = sample(nrow(data), 150)
x_train = data[train_obs, ]
x_holdout = data[-train_obs, ]
outcome_bin = as.numeric(outcome > 22)
y_train = outcome_bin[train_obs]
y_holdout = outcome_bin[-train_obs]
table(y_train, useNA = "ifany")
library(SuperLearner)
listWrappers()
sl_lasso = SuperLearner(Y = y_train, X = x_train, family = binomial(),
                        SL.library = "SL.glmnet")
sl_lasso
sl_lasso$cvRisk[which.min(sl_lasso$cvRisk)]
sl_rf = SuperLearner(Y = y_train, X = x_train, family = binomial(),
                     SL.library = "SL.ranger")
sl = SuperLearner(Y = y_train, X = x_train, family = binomial(),
                  SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"))
pred = predict(sl, x_holdout, onlySL = TRUE)
set.seed(1)

# Don't have timing info for the CV.SuperLearner unfortunately.
# So we need to time it manually.

system.time({
  # This will take about 2x as long as the previous SuperLearner.
  cv_sl = CV.SuperLearner(Y = y_train, X = x_train, family = binomial(),
                          # For a real analysis we would use V = 10.
                          V = 3,
                          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"))
})
summary(cv_sl)











#install.packages("rJava")
#install.packages("bartMachine")
#install.packages("KernelKnn")
#library(rJava)
library(kernlab)
library(KernelKnn)
#library(bartMachine)
data(Boston, package = "MASS")
set.seed(1)
sl_lib = c("SL.xgboost", 
           "SL.randomForest", 
           "SL.glmnet", 
           "SL.nnet", 
           "SL.ksvm",
          "SL.kernelKnn",
          "SL.rpartPrune",
          "SL.lm", 
          "SL.mean")

# Fit XGBoost, RF, Lasso, Neural Net, SVM, BART, K-nearest neighbors, Decision Tree, 
# OLS, and simple mean; create automatic ensemble.
#result = SuperLearner(Y = Boston$medv, X = Boston[, -14], SL.library = sl_lib)

train_obs = sample(nrow(data), 150)
x_train = data[train_obs, ]
x_holdout = data[-train_obs, ]
outcome_bin = outcome
y_train = outcome_bin[train_obs]
y_holdout = outcome_bin[-train_obs]
table(y_train, useNA = "ifany")
library(SuperLearner)
listWrappers()
sl_lasso = SuperLearner(Y = y_train, X = x_train, family = gaussian(),
                        SL.library = sl_lib)
sl_lasso
sl_lasso$cvRisk[which.min(sl_lasso$cvRisk)]
sl_rf = SuperLearner(Y = y_train, X = x_train, family = gaussian(),
                     SL.library = "SL.ranger")
sl = SuperLearner(Y = y_train, X = x_train, family = gaussian(),
                  SL.library = c("SL.glmnet", "SL.ranger"))

pred = predict(sl, x_holdout, onlySL = TRUE)
set.seed(1)
model1 <- lm(y_holdout~pred$pred)
summary(model1)
# Don't have timing info for the CV.SuperLearner unfortunately.
# So we need to time it manually.

system.time({
  # This will take about 2x as long as the previous SuperLearner.
  cv_sl = CV.SuperLearner(Y = y_train, X = x_train, family = binomial(),
                          # For a real analysis we would use V = 10.
                          V = 3,
                          SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"))
})
summary(cv_sl)
