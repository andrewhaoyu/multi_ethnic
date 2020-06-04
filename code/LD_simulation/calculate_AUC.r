#calculate AUC for different ethnic groups
library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01)
#n <- 120000
n.train <- c(100000,15000,15000,15000)
n.test <- c(10000,1500,1500,1500)

#r2 mat represent the r2 matrix for the testing dataset
#column represent the ethnic groups
#row represent different p-value threshold
r2.mat.test <- matrix(0,length(pthres),4)
r2.mat.vad <- matrix(0,length(pthres),4)
for(i in 1:length(eth)){
  #load the phenotype file
  load(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/phenotype.rdata"))
  n <- length(y)
  y.test <- y[(n.train[i]+1):(n.train[i]+n.test[i])]
  y.vad <- y[(n.train[i]+n.test[i]+1):n]
  
  
  
  
  LD <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/LD_clump.clumped")))
  clump.snp <- LD[,3,drop=F]  
  sum.data <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out")))
  colnames(sum.data)[2] <- "SNP"
  prs.all <- left_join(clump.snp,sum.data,by="SNP")
   
  for(k in 1:length(pthres)){
    filedir <- paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/")
    files <- dir(filedir,pattern = paste0("_prs_",k,".profile"))
    prs.score <- rep(0,n)  
    #get the number of
    jdx <- which(prs.all$P<=pthres[k])
    if(length(jdx)>0){
      n.snp.total<-0
      for(j in c(1:22)){
        setwd(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/"))
        filename <- paste0("chr",j,"_prs_",k,".profile")
        idx <- which(prs.all$CHR==j&prs.all$P<=pthres[k])
        n.snp.total=n.snp.total+length(idx)
        if(length(idx)>0){
          prs.temp <- fread(filename)  
          prs.score <- prs.temp$SCORE*2*length(idx)+prs.score
        }
      }
      prs.score = prs.score/(2*n.snp.total)
      #prs.score.mat[,k] = prs.score
      prs.test <- prs.score[(n.train[i]+1):(n.train[i]+n.test[i])]
      prs.vad <- prs.score[(n.train[i]+n.test[i]+1):n]
      #model = lm(y~prs.score)
      model1 <- lm(y.test~prs.test)
      r2.mat.test[k,i] <- summary(model1)$r.square
      model2 <- lm(y.vad~prs.vad)
      r2.mat.vad[k,i] <- summary(model2)$r.square
    }else{
      r2.mat.test[k,i] = 0
      r2.mat.vad[k,i] = 0
      }
    
  }
}
colnames(r2.mat.test) <- colnames(r2.mat.vad) <- eth
rownames(r2.mat.test) <- rownames(r2.mat.vad) <- pthres

r2.mat <- rep(0,length(eth))
#evaluate the best prs performance on the validation
for(k in 1:length(eth)){
  idx <- which.max(r2.mat.test[,k])
  r2.mat[k] <- r2.mat.vad[idx,k]
}
save(r2.mat.test,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/r2.mat.eur.rdata")
write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump.auc.csv")
