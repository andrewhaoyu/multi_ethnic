#Some of the 1KG MEGA data have allele frequencies as 0, need to filter these SNPs out
library(dplyr)
library(data.table)
setwd("/data/zhangh24/KGref_MEGA/GRCh37/")

#load mega SNPs information based on 1KG
mega.snp.infor = fread("/data/zhangh24/multi_ethnic/result/LD_simulation_new/snp_infor_mega+hm3")
mega.snp.infor.select = mega.snp.infor %>% 
  select(SNP, 
         rs_id,
         FREQ_A1_AFR,
         FREQ_A1_AMR,
         FREQ_A1_EUR,
         FREQ_A1_EAS,
         FREQ_A1_SAS)
eth = c("EUR", "AFR", "AMR", "EAS")
i = 1
bim = fread(paste0(eth[i],"/all_chr.bim"))

bim.update = left_join(bim, mega.snp.infor.select, by = c("V2" = "rs_id")) 

var.name = paste0("FREQ_A1_", eth[i])
bim.outlier = bim.update %>% 
  filter(get(var.name) < 0.01 |
          get(var.name) > 0.99)


