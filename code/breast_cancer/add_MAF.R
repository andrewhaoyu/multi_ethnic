#preprocess validation dataset to add MAF infor
setwd("/data/zhangh24/multi_ethnic/data/AABC_data")
load("BC_AFR_overall_valid.rdata")
load("/data/zhangh24/KG.plink/AFR/chr_all_frq.rdata")
freq.infor = freq.infor %>% 
  mutate(chr_pos = paste0(CHR,":",POS))
freq.infor.sub = freq.infor %>% 
  rename(KG.ID = SNP) %>% 
  select(KG.ID,
         chr_pos,MAF)
