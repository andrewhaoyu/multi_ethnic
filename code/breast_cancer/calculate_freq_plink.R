#/data/zhangh24/software/plink2

system(paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/AFR/chr_all --freq --out /data/zhangh24/KG.plink/AFR/chr_all"))
#load 1kg AFR all snps information file
 freq.infor <- fread("/data/zhangh24/KG.plink/AFR/chr_all.frq",header=T)
 freq.infor = freq.infor %>% filter(MAF>=0.005) %>% 
   select(SNP,MAF)
#save(freq.infor,file = "/data/zhangh24/KG.plink/AFR/chr_all_frq.rdata")
library(tidyverse)
freq.infor = freq.infor %>%
  separate(col = SNP,into = c("CHR","POS","A1","A2"),sep = ":",remove =F)
#save(freq.infor,file = "/data/zhangh24/KG.plink/AFR/chr_all_frq.rdata")
#organize 1kg AFR data to have frq with id infor
#load("/data/zhangh24/KG.plink/AFR/chr_all_frq.rdata")

bim <- fread("/data/zhangh24/KG.plink/AFR/chr_all.bim",header=F)
colnames(bim)[1:2] = c("CHR","SNP")
bim = bim %>% select(CHR,SNP)
freq.infor = left_join(freq.infor,bim,by="SNP")

freq.infor = freq.infor %>% 
  rename(CHR = CHR.y) %>% 
  select(SNP,CHR,POS,MAF,A1,A2)
save(freq.infor,file = "/data/zhangh24/KG.plink/AFR/chr_all_frq.rdata")
