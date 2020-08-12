#calculate the AUC for EUR SNPs with EUR coefficients
#with target ethnic coefficients
#with EB coefficients
#calculate AUC for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
#q for three different methods
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
q = as.numeric(args[[4]])

method <- c("eurcoef","tarcoef","eb")

library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01)
#n <- 120000
n.train.vec <- c(15000,45000,80000,100000)
n.train <- n.train.vec[m]
n.test <- (120000-n.train)/2
n.vad <- n.test
n.rep = 100
#r2 mat represent the r2 matrix for the testing dataset
#column represent the ethnic groups
#row represent different p-value threshold
r2.mat.test <- 0
r2.mat.vad <- 0
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
setwd("/data/zhangh24/multi_ethnic/")

#load the phenotype file
y <- as.data.frame(fread(paste0(cur.dir,eth[i],"/phenotypes_rho",l,".phen")))
y <- y[,2+(1:n.rep)]
n <- nrow(y)
y_test_mat <- y[(n.train+1):nrow(y),]
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
load(paste0(cur.dir,"LD.clump.result.rdata"))
load("/data/zhangh24/multi_ethnic/result/LD_simulation/r2.mat.eur.rdata")
#keep the EUR results with sample size at 100,000 (m = 4)
r2.mat <- as.data.frame(LD.result.list[[2]]) %>% 
  filter(eth.vec=="EUR"&
           m_vec==4&
           l_vec==l)
k = which.max(r2.mat$r2.vec)
LD <- as.data.frame(fread(paste0(cur.dir,eth[1],"/LD_clump_rho_",l,"_size_",4,".clumped")))
clump.snp <- LD[,3,drop=F]  
sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_new/",eth[1],"/summary_out_rho_",l,"_size_",m)))  
colnames(sum.data)[2] <- "SNP"
#load EUR best clumped SNP
best.clump.snp <- left_join(clump.snp,sum.data,by="SNP") %>% 
  filter(P<=pthres[k]) %>% 
  select(SNP)

#load target ethnic summary data
sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_new/",eth[i],"/summary_out_rho_",l,"_size_",m)))
prs.all <- inner_join(best.clump.snp,sum.data,by="SNP")

n.snp.total <- 0
prs.score <- rep(0,n)  
for(j in 1:22){
  
#for(k in 1:length(pthres)){
  
  #get the number of
  filename <- paste0(cur.dir,eth[i],"/prs/prs_eursnp_",method[q],"_rho_",l,"_size_",m,"_",j,".profile")
    idx <- which(prs.all$CHR==j) 
    n.snp.total = n.snp.total+length(idx)
    if(length(idx)>0){
      prs.temp <- fread(filename)  
      prs.score <- prs.temp$SCORE*2*length(idx)+prs.score
      }
}
    #prs.score.mat[,k] = prs.score
    prs.test <- prs.score[(n.train+1):(n.train+n.test)]
    prs.vad <- prs.score[(n.train+n.test+1):n]
    #model = lm(y~prs.score)
    r2.test.rep = rep(0,n.rep)
    r2.vad.rep = rep(0,n.rep)
    for(i_rep in 1:n.rep){
      y.test = y_test_mat[1:n.test,i_rep]
      y.vad = y_test_mat[(n.test+1):(nrow(y_test_mat)),i_rep]
      model1 <- lm(y.test~prs.test)
      r2.test.rep[i_rep] <- summary(model1)$r.square
      model2 <- lm(y.vad~prs.vad)
      r2.vad.rep[i_rep]<- summary(model2)$r.square
    }
    r2.mat.test <- mean(r2.test.rep)
    r2.mat.vad <- mean(r2.vad.rep)
  
  
r2 <- mean(c(r2.mat.test,r2.mat.vad))

save(r2,file = paste0(cur.dir,eth[i],"/r2_eursnp_",method[q],"_rho_",l,"_size_",m))
#write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump.auc.csv")
