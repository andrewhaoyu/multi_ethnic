#calculate AUC for different ethnic groups
library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02)
n <- 120000
n.train <- c(100000,15000,15000,15000)
n.test <- c(10000,1500,1500,1500)

#r2 mat represent the r2 matrix for the testing dataset
#column represent the ethnic groups
#row represent different p-value threshold
r2.mat <- matrix(0,length(pthres),4)
for(i in 1:length(eth)){
  #load the phenotype file
  load(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/phenotype.rdata"))
  y.test <- y[(n.train[i]+1):(n.train[i]+n.test[i])]
  n <- length(y)
  prs.score <- rep(0,n)
  for(k in 1:length(pthres)){
    filedir <- paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/")
    files <- dir(filedir,pattern = paste0("_prs_",k,".profile"))
    for(j in c(1:22)){
      setwd(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/"))
      filename <- paste0("chr",j,"_prs_",k,".profile")
      #
      if(filename%in%files){
        prs.temp <- fread(filename)  
        prs.score <- prs.temp$SCORE+prs.score
      }
    }
    prs.test <- prs.score[(n.train[i]+1):(n.train[i]+n.test[i])]
    model1 <- lm(y.test~prs.test)
    r2.mat[k,i] <- summary(model1)$adj.r.squared
    
  }
}
colnames(r2.mat) <- eth
rownames(r2.mat) <- pthres
write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump.auc.csv")
