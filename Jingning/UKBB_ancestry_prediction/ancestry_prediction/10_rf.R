
library(readr)
library(stringr)
pc <- read_tsv("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/plink/chrall.eigenvec")
AFR <- paste0("0:", readLines("/dcs04/nilanjan/data/jzhang2/1000G/1000G_AFR_ID.txt"))
AMR <- paste0("0:", readLines("/dcs04/nilanjan/data/jzhang2/1000G/1000G_AMR_ID.txt"))
EAS <- paste0("0:", readLines("/dcs04/nilanjan/data/jzhang2/1000G/1000G_EAS_ID.txt"))
EUR <- paste0("0:", readLines("/dcs04/nilanjan/data/jzhang2/1000G/1000G_EUR_ID.txt"))
SAS <- paste0("0:", readLines("/dcs04/nilanjan/data/jzhang2/1000G/1000G_SAS_ID.txt"))

## assign true lables for 1000G samples

lab <- paste0(pc$`#FID`,":",pc$IID)
truelabel <- character()
for(i in 1:nrow(pc)){
  if(lab[i] %in% AFR){
    truelabel[i] <- "AFR"
  }else if(lab[i] %in% AMR){
    truelabel[i] <- "AMR"
  }else if(lab[i] %in% EAS){
    truelabel[i] <- "EAS"
  }else if(lab[i] %in% EUR){
    truelabel[i] <- "EUR"
  }else if(lab[i] %in% SAS){
    truelabel[i] <- "SAS"
  }else{
    truelabel[i] <- NA
  }
}

pc <- data.frame(pc)
pc$truelabel <- factor(truelabel)

pc_1000G <- pc[!is.na(pc$truelabel),]

library(caTools)
library(randomForest)

pc_ukbb <- pc[is.na(pc$truelabel),] # those without id are ukbb samples

## some ukbb samples have self-reported labels. Can be used to evaluate accuracy
AFR <- read.table("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/self_reported_AFR.id", header=T)$IID
EAS <- read.table("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/self_reported_EAS.id", header=T)$IID
EUR <- read.table("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/unrelated_whites.id", header=T)$IID
SAS <- read.table("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/self_reported_SAS.id", header=T)$IID

lab <- as.integer(pc_ukbb$IID)
truelabel <- character(length=nrow(pc_ukbb))
m <- which(lab %in% AFR); truelabel[m] <- "AFR"
m <- which(lab %in% EUR); truelabel[m] <- "EUR"
m <- which(lab %in% EAS); truelabel[m] <- "EAS"
m <- which(lab %in% SAS); truelabel[m] <- "SAS"
m <- which(!lab %in% c(AFR,EUR,EAS,SAS)); truelabel[m] <- NA
pc_ukbb$truelabel <- truelabel

########################
## choose number of trees

for(i in c(500,1000,1500,2000,2500,3000)){
set.seed(3)
  classifier_RF = randomForest(x = pc_1000G[c(-1,-2,-ncol(pc_1000G))],
                               y = pc_1000G$truelabel,
                               ntree = i)
saveRDS(classifier_RF, paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/rf_",i,".rds"))

}

#for(i in c(500,1000,1500,2000,2500,3000)){
  i=1500 ## choose the one with maximum accuracy
  set.seed(1)
  y_pred = predict(readRDS(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/rf_",i,".rds")), newdata = pc_ukbb[c(-1,-2,-ncol(pc_ukbb))], type = "prob")

  pred <- character()
  maxp <- numeric()
  for(i in 1:nrow(y_pred)){
    pred[i] <- c("AFR","AMR","EAS","EUR","SAS")[which.max(y_pred[i,])]
    maxp[i] <- max(y_pred[i,])
  }
  table(pred)

  # pred, truelabel
  m <- !is.na(truelabel)
  known_true <- truelabel[m]
  knowm_pred <- pred[m]

  confusion_mtx = table(self_report=known_true, predict=knowm_pred)
  missmatch <- sum(confusion_mtx[1,-1],confusion_mtx[2,-3],confusion_mtx[3,-4],confusion_mtx[4,-5])
  print(missmatch)
#}

#[1] 252
#[1] 255
#[1] 246
#[1] 249
#[1] 252
#[1] 254
#[1] 253
#[1] 255
#[1] 250
#[1] 252


y_pred <- data.frame(y_pred)
y_pred$pred <- pred


pc_ukbb$pred <- pred
#y_pred$true <- pc$truelabel
#m <- which(y_pred$true != y_pred$pred)
#m <- which(classifier_RF$predicted != pc$truelabel)
#y_pred[m,]


table(pred)
 # AFR   AMR   EAS   EUR   SAS
 #9177   787  2020 13902 10858



################################################
## some evaluations (can ignore codes below)

m <- !is.na(truelabel)
table(pred[is.na(truelabel)])

a <- y_pred[is.na(truelabel),]
summary(a[pred[is.na(truelabel)]=="AMR","AMR"])
mean(a[pred[is.na(truelabel)]=="AFR","AFR"]>0.5)
mean(a[pred[is.na(truelabel)]=="AMR","AMR"]>0.5)
mean(a[pred[is.na(truelabel)]=="EAS","EAS"]>0.5)
mean(a[pred[is.na(truelabel)]=="EUR","EUR"]>0.5)
mean(a[pred[is.na(truelabel)]=="SAS","SAS"]>0.5)

a <- y_pred[!is.na(truelabel),]
mean(a[pred[!is.na(truelabel)]=="AFR","AFR"]>0.5)
mean(a[pred[!is.na(truelabel)]=="EAS","EAS"]>0.5)
mean(a[pred[!is.na(truelabel)]=="EUR","EUR"]>0.5)
mean(a[pred[!is.na(truelabel)]=="SAS","SAS"]>0.5)


table(self_report=truelabel[m], predict=pred[m])

maxp <- apply(y_pred,MARGIN=1,max)

final <- data.frame(FID=pc_ukbb$X.FID, IID = pc_ukbb$IID,
                    self_reported_ancestry = pc_ukbb$truelabel, predict_ancestry = pred,
                    voted.p = maxp,
                    flag1.voted.p.0.3 = as.integer(maxp>0.3),
                    flag2.voted.p.0.5 = as.integer(maxp>0.5),
                    flag3.voted.p.0.7 = as.integer(maxp>0.7),
                    flag4.voted.p.0.9 = as.integer(maxp>0.9),
                    p.afr = y_pred[,"AFR"],
                    p.amr = y_pred[,"AMR"],
                    p.eas = y_pred[,"EAS"],
                    p.eur = y_pred[,"EUR"],
                    p.sas = y_pred[,"SAS"]
)

table(final$predict_ancestry)
table(final$predict_ancestry[final$flag1.voted.p.0.3==1])
table(final$predict_ancestry[final$flag2.voted.p.0.5==1])
table(final$predict_ancestry[final$flag3.voted.p.0.7==1])
table(final$predict_ancestry[final$flag4.voted.p.0.9==1])

write_tsv(final,"/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/predicted_ancestry_all.txt")

m <- final$self_reported_ancestry == "EUR"; m[is.na(m)] <- F
write_tsv(final[!m,],"/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/predicted_ancestry_cleaned.txt")


