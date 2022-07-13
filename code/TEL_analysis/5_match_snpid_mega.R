#Goal: match the ID of tel data with MEGA + HM3
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)

eth <- c("EUR","AFR","AMR","EAS","SAS")

setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/TEL_dat/"
#out.dir = paste0("/data/zhangh24/multi_ethnic/result/cleaned/clumping_result/",method,"/",eth[i],"/",trait[l],"/")
load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/snp.infor.match37_38.rdata")
mega.snp = fread("/data/zhangh24/software/PRScsx/1KGLD_MEGA/snpinfo_mult_1kg_hm3")
mega.snp = mega.snp %>% 
  select(SNP) %>% 
  rename(rsid = SNP)
#combine snp_infor file with hm3.snp to get GRCh38 position number
snp.infor = inner_join(snp.infor.match,mega.snp,by=c("rs_id"="rsid")) %>% 
  mutate(chr.pos = paste0(CHR,":",position_GRCh38)) %>% 
  select(rs_id,chr.pos) %>% 
  rename(rsid = rs_id)

#match ID with HM3
for(i in 1:4){
  load(paste0(data.dir,"sum_dat_",eth[i],".rdata")) 
  sum.data = sum_dat %>% 
    mutate(chr.pos = paste0(CHR,":",POS))
  
  #select SNPs exist in HM3 array
  sum.com = inner_join(sum.data,snp.infor,by="chr.pos")
  
  
  sum.com.select = sum.com %>% 
    rename(SNP = rsid,
           A1 = REF,
           A2 = ALT,
           P = PVAL) %>%
    select(CHR,SNP, A1, A2, BETA, P,N,SNP_ID,POS) 
  save(sum.com.select,file = paste0(data.dir,eth[i],"_sum_data_cleaned.rdata"))
  
}

snp.list = list()
for(i in 1:4){
  load(paste0(data.dir,eth[i],"_sum_data_cleaned.rdata"))
  snp.list[[i]] = sum.com.select %>% 
    mutate(V3 = 0) %>% 
    select(CHR,SNP,V3,POS,A1,A2) %>% 
    rename(BP = POS)
  
}
snp.file = rbindlist(snp.list) %>% 
  distinct(SNP,.keep_all=TRUE)
write.table(snp.file,file = paste0(data.dir,"all_eth_mega.bim"),row.names = F,col.names = F,quote=F)
