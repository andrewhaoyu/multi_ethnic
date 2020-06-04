#LD_clumping for different ethnic groups combined with EUR
#LD_clumping for EUR
#load the p-value results
library(data.table)
library(dplyr)
#update the summary results to make it work for plink clumping command
eth <- c("EUR","AFR","AMR","EAS")
#summary.eur <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[1],"/summary.out"),header=T))
summary.eur <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[1],"/summary_MAF.out"),header=T))
colnames(summary.eur)[9] = "peur"
#only keep snpid and p-value for future combination
summary.eur.select = summary.eur %>% 
  select(SNP,peur)
for(i in 2:4){
  #summary <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out"),header=T))
  summary <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary_MAF.out"),header=T))
  #find the shared SNPs between target ethnic group and EUR
  #get the min p-value for between the target ethnic group and EUR for shared snp
  summary.com <- full_join(summary,summary.eur.select,by="SNP")
  summary.com = summary.com %>% 
    mutate(p_update=pmin(P,peur,na.rm = T))
  assoc = summary.com %>%
   # rename(SNP=ID,STAT=T_STAT,BP=POS) %>%
    select(CHR,SNP,BP,A1,TEST,BETA,STAT,p_update) %>% 
    rename(P=p_update)
  write.table(assoc,paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/assoc.out"),col.names = T,row.names = F,quote=F)

}


dim(summary)
head(summary)
pthr = 0.1
r2thr = 0.1
kbpthr = 500
eth <- c("EUR","AFR","AMR","EAS")
code <- rep("c",4)
for(i in 2:4){
  code[i] <- paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/",eth[i],"/chr_all --clump /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/assoc.out --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/LD_clump_two_dim")
  
}
#no need for rerunning in EUR
code <- code[2:4]
write.table(code,file = "/data/zhangh24/multi_ethnic/code/LD_simulation/LD_clumping_two_dim.sh",row.names=F,col.names=F,quote=F)
#the p-value takes the min(AFR,EUR)
#if the SNP only exist


