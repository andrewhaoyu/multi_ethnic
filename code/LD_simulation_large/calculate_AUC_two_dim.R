#calculate AUC for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])

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



n_rep = 10
#for(k in 1:length(pthres)){
#for(i_rep in 1:n.rep){
#for(k in 1:length(pthres)){

summary.eur <- as.data.frame(fread(paste0("./result/LD_simulation_new/",eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep)))  
colnames(summary.eur)[9] = "peur"
colnames(summary.eur)[7] = "beta_eur"
summary.eur.select = summary.eur %>% 
  select(SNP,beta_eur,peur)

#read LD clumped SNPs
LD <- as.data.frame(fread(paste0(cur.dir,eth[i],"/LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,".clumped")))
clump.snp <- LD[,3,drop=F] 
#idx <- which(clump.snp%in%sum.data$SNP==F)
#read the target ethnic group summary level statistics
sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_new/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep))) 
colnames(sum.data)[2] <- "SNP"
#combine the target level summary stat with EUR
summary.com <- left_join(sum.data,summary.eur.select,by="SNP")
#combine the statistics with SNPs after clumping
prs.clump <- left_join(clump.snp,summary.com,by="SNP") 
r2.vec.test <- rep(0,length(pthres)^2)
r2.vec.vad <- rep(0,length(pthres)^2)
pthres1.vec <- rep(0,length(pthres)^2)
pthres2.vec <- rep(0,length(pthres)^2)
temp =1 
for(k1 in 1:length(pthres)){
  for(k2 in 1:length(pthres)){
    
        prs.all <- prs.clump %>% 
          filter((P<=pthres[k1]|
                    peur<=pthres[k2]))
    if(nrow(prs.all)>0){
      n.snp.total <- 0
      prs.score <- rep(0,n)
      for(j in 1:22){
        
        #get the number of
        idx <- which(prs.all$CHR==j)
        n.snp.total = n.snp.total+length(idx)
        if(length(idx)>0){
          
          
          filename <- paste0(cur.dir,eth[i],"/prs/prs_two_dim_",k1,"_",k2,"_rho_",l,"_size_",m,"_chr_",j,"_rep_",i_rep,".profile")
          
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
      #r2.test.rep[i_rep] <- summary(model1)$r.square
      model2 <- lm(y.vad~prs.vad)
      
      r2.vec.test[temp] = summary(model1)$r.square
      r2.vec.vad[temp] = summary(model2)$r.square
      pthres1.vec[temp] = pthres[k1]
      pthres2.vec[temp] = pthres[k2]
      temp = temp + 1
    }else{
      r2.vec.test[temp] = 0
      r2.vec.vad[temp] = 0  
      pthres1.vec[temp] = pthres[k1]
      pthres2.vec[temp] = pthres[k2]
      temp = temp + 1
    }
    
  }
  
}

# }

# }


#}


#names(r2.vec.test) <- names(r2.vec.vad) <- pthres

#evaluate the best prs performance on the validation

idx <- which.max(r2.vec.test)
r2 <- r2.vec.vad[idx]
r2.list <- list(r2,r2.vec.test,r2.vec.vad)
save(r2.list,file = paste0(cur.dir,eth[i],"/r2.list_two_dim_rho_",l,"_size_",m,"_rep_",i_rep))
#write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump.auc.csv")
