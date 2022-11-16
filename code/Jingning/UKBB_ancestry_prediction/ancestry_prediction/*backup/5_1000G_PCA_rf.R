
library(stringr)
pc1000G <- readRDS("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/chr2122/pc1000G_5.pcscores")
AFR <- paste0("0:", readLines("/dcs04/nilanjan/data/jzhang2/1000G/1000G_AFR_ID.txt"))
AMR <- paste0("0:", readLines("/dcs04/nilanjan/data/jzhang2/1000G/1000G_AMR_ID.txt"))
EAS <- paste0("0:", readLines("/dcs04/nilanjan/data/jzhang2/1000G/1000G_EAS_ID.txt"))
EUR <- paste0("0:", readLines("/dcs04/nilanjan/data/jzhang2/1000G/1000G_EUR_ID.txt"))
SAS <- paste0("0:", readLines("/dcs04/nilanjan/data/jzhang2/1000G/1000G_SAS_ID.txt"))

truelabel <- character()
for(i in 1:nrow(pc1000G)){
  if(rownames(pc1000G)[i] %in% AFR){
    truelabel[i] <- "AFR"
  }else if(rownames(pc1000G)[i] %in% AMR){
    truelabel[i] <- "AMR"
  }else if(rownames(pc1000G)[i] %in% EAS){
    truelabel[i] <- "EAS"
  }else if(rownames(pc1000G)[i] %in% EUR){
    truelabel[i] <- "EUR"
  }else if(rownames(pc1000G)[i] %in% SAS){
    truelabel[i] <- "SAS"
  }
}

pc1000G <- data.frame(pc1000G)
pc1000G$truelabel <- factor(truelabel)

# Loading package
library(caTools)
library(randomForest)

set.seed(1)
# Splitting data in train and test data
split <- sample.split(pc1000G, SplitRatio = 0.8)
split

train <- subset(pc1000G, split == "TRUE")
test <- subset(pc1000G, split == "FALSE")

  # Fitting Random Forest to the train dataset
  set.seed(2)  # Setting seed
  classifier_RF = randomForest(x = train[-11],
                               y = train$truelabel,
                               ntree = 500)

  # Predicting the Test set results
  y_pred = predict(classifier_RF, newdata = test[-11])

  # Confusion Matrix
  confusion_mtx = table(test[, 11], y_pred)

#     y_pred
#      AFR AMR EAS EUR SAS
#  AFR 184   0   0   0   0
#  AMR   2  90   0   3   0
#  EAS   0   0 133   0   0
#  EUR   0   1   0 136   0
#  SAS   0   0   0   0 132


y_pred1 = predict(classifier_RF, newdata = test[-11], type = "prob")




set.seed(3)  # Setting seed
  classifier_RF = randomForest(x = pc1000G[-11],
                               y = pc1000G$truelabel,
                               ntree = 500)
saveRDS(classifier_RF, "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/rf.rds")



y_pred1 = predict(classifier_RF, newdata = pc1000G[-11], type = "prob")

pred <- character()
for(i in 1:nrow(y_pred1)){
  pred[i] <- c("AFR","AMR","EAS","EUR","SAS")[which.max(y_pred1[i,])]
}

y_pred1 <- data.frame(y_pred1)
y_pred1$true <- pc1000G$truelabel
y_pred1$pred <- pred
m <- which(y_pred1$true != y_pred1$pred)

m <- which(classifier_RF$predicted != pc1000G$truelabel)
y_pred1[m,]

