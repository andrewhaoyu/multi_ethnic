#calculate AUC for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
k = as.numeric(args[[4]])

library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
#n <- 120000
n.train.vec <- c(15000,45000,80000,100000)
n.train <- n.train.vec[m]
n.test <- (120000-n.train)/2
n.vad <- n.test
n.rep = 10
#r2 mat represent the r2 matrix for the testing dataset
#column represent the ethnic groups
#row represent different p-value threshold
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
setwd("/data/zhangh24/multi_ethnic/")

  #load the phenotype file
  y <- as.data.frame(fread(paste0(cur.dir,eth[i],"/phenotypes_rho",l,".phen")))
  y <- y[,2+(1:n.rep)]
  n <- nrow(y)
  y_test_mat <- y[(n.train+1):nrow(y),]
  #y.test <- y[(n.train[i]+1):(n.train[i]+n.test[i])]
  # y.vad <- y[(n.train[i]+n.test[i]+1):n]
  #
  
  
  
  LD <- as.data.frame(fread(paste0(cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,".clumped")))
  clump.snp <- LD[,3,drop=F]  
  sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_new/",eth[i],"/summary_out_rho_",l,"_size_",m)))  
  colnames(sum.data)[2] <- "SNP"
  
  #for(k in 1:length(pthres)){
  r2.vec.test <- rep(0,i_rep)
  r2.vec.vad <- rep(0,i_rep)
  for(i_rep in 1:n.rep){
  #for(k in 1:length(pthres)){
    n.snp.total <- 0
    prs.score <- rep(0,n)  
    prs.all <- left_join(clump.snp,sum.data,by="SNP")  %>% 
    filter(P<=pthres[k])
    if(nrow(prs.all)>0){
      for(j in 1:22){
        
        #get the number of
        idx <- which(prs.all$CHR==j)
        n.snp.total = n.snp.total+length(idx)
        if(length(idx)>0){
          
          
          filename <- paste0(cur.dir,eth[i],"/prs/prs_rho_",l,"_size_",m,"_chr_",j,"_rep_",i_rep,".profile")
          
          prs.temp <- fread(filename)  
          prs.score <- prs.temp$SCORE*2*length(idx)+prs.score
        }
      }
      #prs.score.mat[,k] = prs.score
      prs.test <- prs.score[(n.train+1):(n.train+n.test)]
      prs.vad <- prs.score[(n.train+n.test+1):n]
      #model = lm(y~prs.score)
      y.test = y_test_mat[1:n.test,i_rep]
      y.vad = y_test_mat[(n.test+1):(nrow(y_test_mat)),i_rep]
      model1 <- lm(y.test~prs.test)
      r2.test.rep[i_rep] <- summary(model1)$r.square
      model2 <- lm(y.vad~prs.vad)
      r2.mat.test[i_rep,k] = summary(model1)$r.square
      r2.mat.vad[i_rep,k] = summary(model2)$r.square
    }
    r2.mat.test[i_rep,k] = 0
    r2.mat.vad[i_rep,k] = 0
   # }
    
  }
      
   
    #}


names(r2.mat.test) <- names(r2.mat.vad) <- pthres

#evaluate the best prs performance on the validation

  idx <- which.max(r2.mat.test)
  r2 <- r2.mat.vad[idx]
r2.list <- list(r2,r2.mat.test,r2.mat.vad)
save(r2.list,file = paste0(cur.dir,eth[i],"/r2.list_rho_",l,"_size_",m))
#write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump.auc.csv")
