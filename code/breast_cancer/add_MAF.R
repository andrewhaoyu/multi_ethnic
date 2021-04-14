#preprocess validation dataset to add MAF infor
setwd("/data/zhangh24/multi_ethnic/data/AABC_data")
load("/data/zhangh24/KG.plink/AFR/chr_all_frq.rdata")
#library(dplpyr)
library(tidverse)
library(data.table)
freq.infor = freq.infor %>% 
  unite(chr_pos_A1_A2,CHR,POS,A1,A2,sep=":",remove=F) %>% 
  unite(chr_pos_A2_A1,CHR,POS,A2,A1,sep=":",remove=F)
freq.infor.sub = freq.infor %>% 
  rename(KG.ID = SNP,
         KG.A1 = A1,
         KG.A2 = A2) %>%
  select(KG.ID,
         chr_pos_A1_A2,
         chr_pos_A2_A1)


load("BC_AFR_overall_valid.rdata")
#library(fuzzyjoin)
library(tidyverse)
# freq.infor.sub = freq.infor.sub %>% 
#   unite(chr_pos_unit,chr_pos_A1_A2,chr_pos_A2_A1,remove=F)
# sum.data = sum.data %>% 
#   mutate(chr_pos = paste0(CHR,":",POS))
sum.data.update = sum.data %>% 
  mutate(chr_pos = paste0(CHR,":",POS))
sum.data.update = sum.data %>% 
  left_join(freq.infor.sub,by=c("ID" = "chr_pos_A1_A2")) %>% 
  select(c(colnames(sum.data),"KG.ID")) %>% 
  left_join(freq.infor.sub,by=c("ID" = "chr_pos_A2_A1")) %>% 
  mutate(KG.ID = ifelse(is.na(KG.ID.x),KG.ID.y,KG.ID.x)) %>% 
  select(c(colnames(sum.data),"KG.ID"))

save(sum.data.update,file = "/data/zhangh24/multi_ethnic/data/AABC_data/BC_AFR_overall_valid_KGID.rdata")





load("BC_AFR_overall_train.rdata")
sum.data = sum.data %>% 
  mutate(chr_pos = paste0(CHR,":",POS))
sum.data.update = sum.data %>% 
  left_join(freq.infor.sub,by=c("ID" = "chr_pos_A1_A2")) %>% 
  select(c(colnames(sum.data),"KG.ID")) %>% 
  left_join(freq.infor.sub,by=c("ID" = "chr_pos_A2_A1")) %>% 
  mutate(KG.ID = ifelse(is.na(KG.ID.x),KG.ID.y,KG.ID.x)) %>% 
  select(c(colnames(sum.data),"KG.ID"))

save(sum.data.update,file = "/data/zhangh24/multi_ethnic/data/AABC_data/BC_AFR_overall_train_KGID.rdata")







#subset the breast cancer data to mega list
#subset the breast cancer data to mega list + best eur prs
#bim <- as.data.frame(fread("/data/zhangh24/KGref_MEGA/GRCh37/AFR/all_chr.bim"))
#prs.snp <- read.csv("/data/zhangh24/multi_ethnic/data/breast_cancer_330SNPs.csv",header=T)

#idx <- which(prs.snp$variant%in%bim$V2==F)






# 
# mega.list <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_new/mega-hm3-rsid.txt"),header  =F))
# 
# load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/snp.infor.match37_38.rdata")
