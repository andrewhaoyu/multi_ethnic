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
np = nrow(expand.grid(pthres,pthres))
r2.mat <- cbind(expand.grid(pthres,pthres),0,0,0)
r2.mat.test <- cbind(expand.grid(pthres,pthres),0,0,0)
r2.mat.vad <- cbind(expand.grid(pthres,pthres),0,0,0)

# summary.eur <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[1],"/summary.out"),header=T))
# colnames(summary.eur)[9] = "peur"
# colnames(summary.eur)[7] = "beta_eur"
# summary.eur.select = summary.eur %>% 
#   select(SNP,beta_eur,peur)
n.rep = 100
for(i in 2:length(eth)){
  #load the phenotype file
  load(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/phenotype.rdata"))
  n <- length(y)
  # y.test <- y[(n.train[i]+1):(n.train[i]+n.test[i])]
  # y.vad <- y[(n.train[i]+n.test[i]+1):n]
  load(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/phenotype_test_mat.rdata"))
  load(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/LD_clump_two_dim.clumped_all_infor"))
  #k for the p-value threshold on the target population
  #l for the p-value threshold on the eur population
  
  for(k in 1:length(pthres)){
    for(l in 1:length(pthres)){
      prs.score <- rep(0,n)
      prs.file <- prs.all %>% filter(P<=pthres[k]|
                                       peur<=pthres[l])%>% 
        select(SNP,A1,BETA,CHR)
      
      filedir <- paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/")
      files <- dir(filedir,pattern = paste0("_prs_",k,".profile"))
      #get the number of
      kdx <- nrow(prs.file)
      if(length(kdx)>0){
        n.snp.total<-0
        for(j in c(1:22)){
          setwd(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/"))
          filename <- paste0("chr",j,"_prs_eb_",k,"_",l,".profile")
          idx <- which(prs.file$CHR==j)
          n.snp.total=n.snp.total+length(idx)
          if(length(idx)>0){
            prs.temp <- fread(filename)  
            prs.score <- prs.temp$SCORE*2*length(idx)+prs.score
          }
        }
        prs.score = prs.score/(2*n.snp.total)
        prs.test <- prs.score[(n.train[i]+1):(n.train[i]+n.test[i])]
        prs.vad <- prs.score[(n.train[i]+n.test[i]+1):n]
        r2.test.rep = rep(0,n.rep)
        r2.vad.rep = rep(0,n.rep)
        for(m in 1:n.rep){
          y.test = y_test_mat[1:n.test[i],m]
          y.vad = y_test_mat[(n.test[i]+1):(nrow(y_test_mat)),m]
          model1 <- lm(y.test~prs.test)
          r2.test.rep[m] <- summary(model1)$r.square
          model2 <- lm(y.vad~prs.vad)
          r2.vad.rep[m]<- summary(model2)$r.square
        }
        jdx = which(r2.mat.test[,1]==pthres[k]&r2.mat.test[,2]==pthres[l])
        r2.mat.test[jdx,i+1] <- mean(r2.test.rep)
        r2.mat.vad[jdx,i+1] <- mean(r2.vad.rep)
        
        # prs.test <- prs.score[(n.train[i]+1):(n.train[i]+n.test[i])]
        # prs.vad <- prs.score[(n.train[i]+n.test[i]+1):n]
        # model1 <- lm(y.test~prs.test)
        # jdx = which(r2.mat.test[,1]==pthres[k]&r2.mat.test[,2]==pthres[l])
        # r2.mat.test[jdx,i+1] <- summary(model1)$r.square
        # model2 <- lm(y.vad~prs.vad)
        # jdx = which(r2.mat.vad[,1]==pthres[k]&r2.mat.vad[,2]==pthres[l])
        # r2.mat.vad[jdx,i+1] <- summary(model1)$r.square
      }else{
        r2.mat.test[jdx,i+1]=0
        r2.mat.vad[jdx,i+1] = 0
      }
      
    }
    
  }
}
colnames(r2.mat.test)[3:5] <- colnames(r2.mat.vad)[3:5] <- eth[2:4]
r2.max <- rep(0,3)
for(i in 1:3){
  jdx <- which.max(r2.mat.test[,i+2])
  r2.max[i] = max(r2.mat.vad[jdx,i+2])  
}


write.csv(r2.max,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/r2.max_eb.csv")
