#use plink2 to calculate prs
#load LD_clump_file
library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01)
n.snp.mat <- matrix(0,length(pthres),4)
for(i in 1:length(eth)){
  LD <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/LD_clump.clumped")))
  clump.snp <- LD[,3,drop=F]  
  sum.data <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out")))
  colnames(sum.data)[2] <- "SNP"
  prs.all <- left_join(clump.snp,sum.data,by="SNP") 
  for(j in 1:length(pthres)){
    prs.file <- prs.all %>% filter(P<=pthres[j]) %>% 
      select(SNP,A1,BETA)
    n.snp.mat[j,i] <- nrow(prs.file)
    write.table(prs.file,file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_file_pvalue_",j),col.names = T,row.names = F,quote=F)
    }
  
}
write.csv(n.snp.mat,file =paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/n_snp_mat.csv"))
# i <- 4
# n.train <- 15000
# load("/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/phenotype.rdata")
# load("/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/causal_genotype.rdata")
# sum.data <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out")))
# model <- lm(y[1:n.train[i]]~genotype[jdx,1:n.train[i]])
# 
# summary(model)
# model <- lm(y~genotype[jdx,])
# idx <- which(sum.data$SNP=="7:32777664:G:A")
# jdx <- which(row.names(genotype)=="7:32777664:G:A")
# prs <- genotype[jdx,]*-0.228463
# n.train <- c(100000,15000,15000,15000)
# n.test <- c(10000,1500,1500,1500)
# 
# model2 <- lm(y~prs)
# summary(model2)
# 
# y.test <- y[(n.train[i]+1):(n.train[i]+2*n.test[i])]
# prs.test <- prs[(n.train[i]+1):(n.train[i]+2*n.test[i])]
# 
# model2 <- lm(y~prs-1)
# summary(model2)
# model3 <- lm(y~prs)
# summary(model3)
# 
# sum(genotype[jdx,])/(2*ncol(genotype))
# p <- rep(0,nrow(genotype))
# for(k in 1:length(p)){
#   print(k)
#   model <- lm(y[1:n.train]~genotype[k,1:n.train])
#   p[k] <- coefficients(summary(model))[2,4]
# }

#use plink score funciton to get prs
code <- rep("c",10000)
temp <- 1
for(i in 1:length(eth)){
  for(j in 1:22){
    for(k in 1:length(pthres)){
      temp.code <- paste0("/data/zhangh24/software/plink2 --score /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_file_pvalue_",k," no-sum no-mean-imputation --bfile /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/chr",j,"_prs_",k)
      code[temp] <- temp.code
      temp <- temp+1
    }
  }
  
}
code <- code[1:(temp-1)]
write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/calculate_prs.sh"),col.names = F,row.names = F,quote=F)







