library(dplyr)
library(data.table)
setwd("/data/zhangh24/multi_ethnic/result/TEL/prscsx_mega")
eth = c("EUR", "AFR", "AMR", "EAS")
phi = c("1e+00", "1e-02", "1e-04", "1e-06","1e-08")
data.dir = "/data/zhangh24/multi_ethnic/data/TEL_dat/"
for(i in 1:4){
  load(paste0(data.dir,eth[i],"_sum_data_cleaned_mega.rdata"))
  snp.id = sum.com.select %>% 
    select(SNP, SNP_ID)
  for(v in 1:5){
    snp.list = list()
    for(j in 1:22){
      print(j)
      data = fread(paste0("sum_",eth[i],"_pst_eff_a1_b0.5_phi",phi[v],"_chr",j,".txt"))
      colnames(data) = c("CHR", "rsid", "POS", "effect_allele", "non_effect_allele", "post_effect")
      data.com = left_join(data, snp.id, by = c("rsid"="SNP"))
      snp.list[[j]] = data.com        
    }
    cleaned.result = rbindlist(snp.list)
    save(cleaned.result, file = paste0("./prs_csx_mega_result_cleaned/sum_",eth[i],"_pst_eff_a1_b0.5_phi",phi[v],".rdata"))
  }
}



