#LD_clumping for different ethnic groups
#LD_clumping for EUR
#load the p-value results
library(data.table)
library(dplyr)
#update the summary results to make it work for plink clumping command
eth <- c("EUR","AFR","AMR","EAS")
for(i in 1:4){
  summary <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out"),header=T))
  colnames(summary)[1] <- "CHR"
  assoc = summary %>% 
    rename(SNP=ID,STAT=T_STAT,BP=POS) %>% 
    select(CHR,SNP,BP,A1,TEST,BETA,STAT,P)
  write.table(assoc,paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/assoc.out"),col.names = T,row.names = F,quote=F)
  
}


dim(summary)
head(summary)
pthr = 0.01
r2thr = 0.1
kbpthr = 500
eth <- c("EUR","AFR","AMR","EAS")
code <- rep("c",4)
for(i in 1:4){
  code[i] <- paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/",eth[i],"/chr_all --clump /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/assoc.out --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/LD_clump")
  
}
write.table(code,file = "/data/zhangh24/multi_ethnic/code/LD_simulation/LD_clumping.sh",row.names=F,col.names=F,quote=F)
#the p-value takes the min(AFR,EUR)
#if the SNP only exist