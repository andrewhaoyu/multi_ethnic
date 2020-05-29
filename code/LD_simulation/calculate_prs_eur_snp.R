#use plink2 to calculate for EUR best PRS, use EUR cofficients or AFR cofficients
#load LD_clump_file

library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
n.snp.mat <- matrix(0,length(pthres),4)
#load EUR r2 result
load("/data/zhangh24/multi_ethnic/result/LD_simulation/r2.mat.eur.rdata")
#get the best performance eur prs p-value threshold
jdx = which.max(r2.mat[,1])
i= 1
LD <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/LD_clump.clumped")))
clump.snp <- LD[,3,drop=F]  
sum.data <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out")))
colnames(sum.data)[2] <- "SNP"
prs.all <- left_join(clump.snp,sum.data,by="SNP") 
n.snp.mat <- matrix(0,1,4)
#use EUR regression coefficients
for(i in 2:length(eth)){
  
    prs.file <- prs.all %>% filter(P<=pthres[jdx]) %>% 
      select(SNP,A1,BETA)
    #n.snp.mat[1,i] <- nrow(prs.file)
    write.table(prs.file,file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_eursnp_eurcoef"),col.names = T,row.names = F,quote=F)
  
  
}
#use target population regresssion coefficients
#get the number of SNPs in different population
n.snp.total = rep(0,4)
n.snp.total[1] = nrow(prs.file)
for(i in 2:length(eth)){
  prs.file <- prs.all %>% filter(P<=pthres[jdx]) %>% 
    select(SNP)
  sum.data <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out")))
  #combine the prs with the target population coefficients
  prs.file.com <- left_join(prs.file,sum.data,by="SNP")
  prs.file <- prs.file.com %>% select(SNP,A1,BETA)
  prs.file <- prs.file[complete.cases(prs.file),]
  n.snp.total[i] = nrow(prs.file)
  #n.snp.mat[j,i] <- nrow(prs.file)
  write.table(prs.file,file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_eursnp_tarcoef"),col.names = T,row.names = F,quote=F)
  
  
}
#use plink score funciton to get prs
code <- rep("c",10000)
temp <- 1
for(i in 2:length(eth)){
  for(j in 1:22){
      temp.code <- paste0("/data/zhangh24/software/plink2 --score /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_eursnp_eurcoef no-sum no-mean-imputation --bfile /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/chr",j,"_prs_eursnp_eurcoef")
      code[temp] <- temp.code
      temp <- temp+1
    
  }
  
}


for(i in 2:length(eth)){
  for(j in 1:22){
    temp.code <- paste0("/data/zhangh24/software/plink2 --score /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_eursnp_tarcoef no-sum no-mean-imputation --bfile /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/chr",j,"_prs_eursnp_tarcoef")
    code[temp] <- temp.code
    temp <- temp+1
    
  }
  
}

code <- code[1:(temp-1)]
write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/calculate_prs_eursnp.sh"),col.names = F,row.names = F,quote=F)







