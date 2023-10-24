#load the summary statistics

library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")

#load summary level statistics for overall
#Haoyupdate added the meta-analysis between iCOGs and OncoArray
standard_result <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF_Haoyuupdate.txt"))
colnames(standard_result)[1] <- "MarkerName"

bcac_result = fread("/data/zhangh24/ldsc/bcac_result.txt")
bcac_result = bcac_result[,-c(1:2)]
write.table(bcac_result,file = "/data/zhangh24/multi_ethnic/result/stat_gene_course/data/overall_bc",
            row.names = F, col.names = T, quote = F)


load(paste0("./genetic_correlation/result/BCAC_subtypes_result_082119.Rdata"))
head(BCAC_subtypes_result)
load("./genetic_correlation/result/hapmap3list.Rdata")
shared.data = shared.data %>% select(SNP,SNP.ICOGS)
library(dplyr)
update_data = left_join(shared.data,BCAC_subtypes_result, by = c("SNP.ICOGS"="rs_id"))
lua_sum <- update_data %>% 
  select(SNP,
         CHR,
         position,
         exp_freq_a1,
         references_allele,
         effect_allele,
         log_or_Luminial_A,
         var_Luminial_A,
         z_stat_Luminial_A,
         effect_sample_size_Luminial_A,
         p_value_Luminial_A) %>% 
  mutate(snpid = SNP,
         bp = position,
         A1 = effect_allele,
         A2 = references_allele,
         Z = z_stat_Luminial_A,
         P = p_value_Luminial_A,
         info = 1,
         MAF = ifelse(exp_freq_a1<=0.5,exp_freq_a1,1-exp_freq_a1),
         N = effect_sample_size_Luminial_A) %>% 
  select(snpid, CHR, bp, A2, A1, Z, P, info, MAF, N)

write.table(lua_sum,file = "/data/zhangh24/multi_ethnic/result/stat_gene_course/data/lua_bc",
            row.names = F, col.names = T, quote = F)


tn_sum <- update_data %>% 
  select(SNP,
         CHR,
         position,
         exp_freq_a1,
         references_allele,
         effect_allele,
         log_or_Triple_Negative,
         var_Triple_Negative,
         z_stat_Triple_Negative,
         effect_sample_size_Triple_Negative,
         p_value_Triple_Negative) %>% 
  mutate(snpid = SNP,
         bp = position,
         A1 = effect_allele,
         A2 = references_allele,
         Z = z_stat_Triple_Negative,
         P = p_value_Triple_Negative,
         info = 1,
         MAF = ifelse(exp_freq_a1<=0.5,exp_freq_a1,1-exp_freq_a1),
         N = effect_sample_size_Triple_Negative) %>% 
  select(snpid, CHR, bp, A2, A1, Z, P, info, MAF, N)

write.table(tn_sum,file = "/data/zhangh24/multi_ethnic/result/stat_gene_course/data/tn_bc",
            row.names = F, col.names = T, quote = F)
