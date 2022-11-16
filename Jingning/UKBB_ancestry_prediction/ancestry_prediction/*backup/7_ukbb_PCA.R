

suppressMessages(library('plink2R'))

genos = read_plink(paste("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_final/chr2122"), impute="avg")
#a <- colnames(genos$bed)
#b <- rownames(loadings)
#identical(a,b)


## memory cannot map
#genos$bed <- scale(genos$bed)

#ukb_sd <- apply(genos$bed, MARGIN=2, sd)
#ukb_mean <- apply(genos$bed, MARGIN=2, mean)
#save(ukb_sd, ukb_mean, file="/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/hdpca/ukbbPC_original/geno.sdmean.RData")

load("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/hdpca/geno1000G.sdmean.RData")
for(i in 1:ncol(genos$bed)){
  genos$bed[,i] <- (genos$bed[,i] - genomes_mean[i]) / genomes_sd[i]
  if(i %% 1000 == 0){
    print(i)
  }
}
saveRDS(genos, "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/hdpca/ukbbPC_original/std_geno.rds")

loadings <- readRDS( "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/chr2122/pc1000G_5.loadings")

#m <- which(ukb_sd==0)
#genos$bed <- genos$bed[,-m]
#loadings <- loadings[-m,]

## memory cannot map
#ukbb_pc <- genos$bed %*% loadings

ncol(genos$bed)/1000
for(i in 1:84){
  print(i)
  if( i ==1 ){
    res <- genos$bed[,((i-1)*1000+1):(i*1000)] %*% loadings[((i-1)*1000+1):(i*1000),]
  }else if(i > 1 & i != 84){
    res <- res + genos$bed[,((i-1)*1000+1):(i*1000)] %*% loadings[((i-1)*1000+1):(i*1000),]
  }else{
    res <- res + genos$bed[,((i-1)*1000+1):ncol(genos$bed)] %*% loadings[((i-1)*1000+1):ncol(genos$bed),]
  }
}
ukbb_pc <- res

saveRDS(ukbb_pc, "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/hdpca/ukbbPC_original/ukbb_pc.rds")


library(hdpca)
load("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/hdpca/parameters.RData")
ukbb_pc.adj <- pc_adjust(train.eval,p,n,ukbb_pc,method="d.gsp",n.spikes=m)

library(randomForest)
classifier_RF <- readRDS("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/rf.rds")

y_pred = predict(classifier_RF, newdata = ukbb_pc, type = "prob")

pred <- character()
prob <- numeric()
for(i in 1:nrow(y_pred)){
  pred[i] <- c("AFR","AMR","EAS","EUR","SAS")[which.max(y_pred[i,])]
  prob[i] <- max(y_pred[i,])
}

y_pred <- data.frame(y_pred)
y_pred$maxprob <- prob
y_pred$predicted <- pred
#saveRDS(y_pred, "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/predicted_ancestry/y_pred.rds")
y_pred$IID <- unlist(lapply(str_split(rownames(y_pred), ":"), FUN=function(x){x[2]}))
table(y_pred$predicted)

EAS <- fread2("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/self_reported_EAS.id")
a <- y_pred[y_pred$IID %in% EAS$IID,]
head(a)

nrow(y_pred)
mean(y_pred$maxprob>0.5)





