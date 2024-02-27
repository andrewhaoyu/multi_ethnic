
library(data.table)
library(dplyr)
snp_infor = readRDS("/data/BB_Bioinformatics/ProjectData/SNP_infor_GRCh37_38/SNP_GRCh37_38_match_update.rds")
snp_infor_sub = snp_infor %>% select(rsid,pos37)
library(dplyr)
trait_vec = c("height","bmi")
eth_vec = c("EUR","AFR")
method_vec = c("LDpred2","weighted_LDpred2")
setwd("/data/zhangh24/multi_ethnic/result/PGS_LDpred2score")
files = dir(pattern = ".txt.gz",path = "/data/zhangh24/multi_ethnic/result/PGS_LDpred2score")
for(i in 1:2){
  for(l in 1:2){
    for(m in 1:2){
      trait = trait_vec[i]
      eth = eth_vec[l]
      method = method_vec[m]
      filename = paste0(trait,"_",eth,"_",method,".txt.gz")
      if(filename %in% files){
        data = read.table(gzfile(filename),header = T)
        data_update = left_join(data,snp_infor_sub,by = c("rsID"="rsid"))
        sum(is.na(data_update$pos37))
        data_update_select = data_update %>% 
          select(rsID,chr_name,pos37,effect_allele,
                 other_allele,effect_weight) %>% 
          rename(chr_position = pos37)
       write.table(data_update_select,
                   file = gzfile(paste0("./PGS_LDpred2score_update/",filename)),
       row.names = F, col.names = T, quote = F)
      }
      
    }
  }
}
