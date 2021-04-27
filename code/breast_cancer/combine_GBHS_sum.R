
setwd("/data/zhangh24/multi_ethnic/data/GBHS_sum/")

pheno <- c("All_cases_Nov_1","ER-_cases_Nov_1","ER+_cases_Nov_1")
outpheno <- c("overall","erpos","erneg")
library(data.table)
library(dplyr)

for(i1 in 1:3){
  data.list = list()
    for(j in 1:22){
    data = fread(paste0(pheno[i1],"/","summary_All_sub_type_add_info_",j,".csv"),header=T)  
    data.select = data %>% 
      rename(CHR = chromosome,
             POS = position,
             Ref_allele = alleleA,
             Eff_allele = alleleB,
             MAF = controls_maf,
             BETA = frequentist_add_beta_1,
             SE = frequentist_add_se_1,
             P = frequentist_add_pvalue) %>% 
      select(CHR,POS,Ref_allele,
             Eff_allele,MAF,BETA,
             SE,P)
    data.list[[j]] = data.select
  }
  sum.data = rbindlist(data.list)
  save(sum.data, file = paste0(outpheno[i1],"_sum.rdata"))
  
}
