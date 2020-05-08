#use plink2 to calculate prs
#load LD_clump_file
library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02)

for(i in 1:length(eth)){
  LD <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/LD_clump.clumped")))
  clump.snp <- LD[,3,drop=F]  
  sum.data <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out")))
  colnames(sum.data)[3] <- "SNP"
  prs.all <- left_join(clump.snp,sum.data,by="SNP") 
  for(j in 1:length(pthres)){
    prs.file <- prs.all %>% filter(P<=pthres[j]) %>% 
      select(SNP,A1,BETA)
    write.table(prs.file,file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_file_pvalue_",j),col.names = T,row.names = F,quote=F)
    }
  
}



#use plink score funciton to get prs
code <- rep("c",10000)
temp <- 1
for(i in 1:length(eth)){
  for(j in 1:22){
    for(k in 1:length(pthres)){
      temp.code <- paste0("/data/zhangh24/software/plink2 --score /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_file_pvalue_",k," no-sum no-mean-imputation --gen /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".plink.tag.gen --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id --sample /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/sample.txt --out /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/chr",j,"_prs_",k)
      code[temp] <- temp.code
      temp <- temp+1
    }
  }
  
}
code <- code[1:(temp-1)]
write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/calculate_prs.sh"),col.names = F,row.names = F,quote=F)







