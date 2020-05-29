#calculate AUC for different ethnic groups use EUR prs
library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS")
#n <- 120000
n.train <- c(100000,15000,15000,15000)
n.test <- c(10000,1500,1500,1500)

#r2 mat represent the r2 matrix for the testing dataset
#column represent the ethnic groups
#row represent different p-value threshold
r2.mat <- matrix(0,1,4)
i = 1
LD <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/LD_clump.clumped")))
clump.snp <- LD[,3,drop=F]  
sum.data <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out")))
colnames(sum.data)[2] <- "SNP"
prs.all <- left_join(clump.snp,sum.data,by="SNP")
#get the best eur prs SNPs
colnames(sum.data)[2] <- "SNP"
prs.all <- left_join(clump.snp,sum.data,by="SNP") 
#load EUR r2 result
load("/data/zhangh24/multi_ethnic/result/LD_simulation/r2.mat.eur.rdata")
#get the best performance eur prs p-value threshold
jdx = which.max(r2.mat[,1])
prs.eur.snp <- prs.all %>% filter(P<=pthres[jdx]) 

for(i in 2:length(eth)){
  #load the phenotype file
  load(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/phenotype.rdata"))
  y.test <- y[(n.train[i]+1):(n.train[i]+n.test[i])]
  n <- length(y)
  
  
  
  prs.score <- rep(0,n)  
   
    
      n.snp.total<-0
      for(j in c(1:22)){
        setwd(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/"))
        filename <- paste0("chr",j,"_prs_eursnp_eurcoef.profile")
        idx <- which(prs.eur.snp$CHR==j)
       n.snp.total=n.snp.total+length(idx)
        if(length(idx)>0){
          prs.temp <- fread(filename)  
          prs.score <- prs.temp$SCORE*2*length(idx)+prs.score
        }
      }
      prs.score = prs.score/(2*n.snp.total)
      #prs.score.mat[,k] = prs.score
      prs.test <- prs.score[(n.train[i]+1):(n.train[i]+n.test[i])]
      
      model1 <- lm(y.test~prs.test)
      r2.mat[1,i] <- summary(model1)$r.square
    }
    

colnames(r2.mat) <- eth
#rownames(r2.mat) <- pthres
#save(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/r2.mat.eur.rdata")
write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump_eurprseurcoef.csv")



for(i in 2:length(eth)){
  #load the phenotype file
  load(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/phenotype.rdata"))
  y.test <- y[(n.train[i]+1):(n.train[i]+n.test[i])]
  n <- length(y)
  
  
  
  prs.score <- rep(0,n)  
  
  
  n.snp.total<-0
  for(j in c(1:22)){
    setwd(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/"))
    filename <- paste0("chr",j,"_prs_eursnp_tarcoef.profile")
    idx <- which(prs.eur.snp$CHR==j)
    n.snp.total=n.snp.total+length(idx)
    if(length(idx)>0){
      prs.temp <- fread(filename)  
      prs.score <- prs.temp$SCORE*2*length(idx)+prs.score
    }
  }
  prs.score = prs.score/(2*n.snp.total)
  #prs.score.mat[,k] = prs.score
  prs.test <- prs.score[(n.train[i]+1):(n.train[i]+n.test[i])]
  
  model1 <- lm(y.test~prs.test)
  r2.mat[1,i] <- summary(model1)$r.square
}


colnames(r2.mat) <- eth
