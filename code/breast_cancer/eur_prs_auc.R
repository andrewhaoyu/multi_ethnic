setwd("/data/zhangh24/multi_ethnic/data/AABC_data")
load("BC_AFR_overall_train.rdata")
prs.snp <- read.csv("/data/zhangh24/multi_ethnic/data/breast_cancer_330SNPs.csv",header=T)
var_name = data.frame(ID = gsub("_",":",prs.snp[,1]),rsID = prs.snp[,2])
# idx <- which(var_name%in%sum.data[,1]==F)
# 
# jdx <- which(sum.data$CHR==17&
#                sum.data$POS==7571752)
library(dplyr)
best.eur.snp = inner_join(var_name,sum.data,
                         by="ID",)

#prepare data for summary AUC
prs.model.file = best.eur.snp %>% 
  mutate(Effect_Allele = ALT,
         BETA = beta.ALT) %>% 
  select(rsID,Effect_Allele,BETA)

write.table(prs.model.file,
            file = "/data/zhangh24/multi_ethnic/result/breast_cancer/best_eur_prsmodelfile",
            row.names = F,
            col.names = F,
            quote=F)
load("BC_AFR_overall_valid.rdata")

best.eur.snp.select = best.eur.snp %>% 
  select(ID,rsID,ALT) %>% 
  rename(A1 = ALT)

gwas.summary.data.test = left_join(best.eur.snp.select,
                                   sum.data,by="ID")
all.equal(gwas.summary.data.test$A1,
          gwas.summary.data.test$ALT)
library(data.table)
#load 1kg AFR all snps information file
# freq.infor <- fread("/data/zhangh24/KG.plink/AFR/chr_all.frq",header=T)
# freq.infor = freq.infor %>% filter(MAF>=0.005) %>% 
#   select(SNP,MAF)
#save(freq.infor,file = "/data/zhangh24/KG.plink/AFR/chr_all_frq.rdata")
#library(tidyverse)
# freq.infor = freq.infor %>% 
#   separate(col = SNP,into = c("CHR","POS","A1","A2"),sep = ":",remove =F)
#save(freq.infor,file = "/data/zhangh24/KG.plink/AFR/chr_all_frq.rdata")
load("/data/zhangh24/KG.plink/AFR/chr_all_frq.rdata")
freq.infor = freq.infor %>% 
  select(POS,MAF) %>% 
  mutate(POS = as.numeric(POS))

gwas.summary.data.test = left_join(gwas.summary.data.test,freq.infor,by="POS")


gwas.summary.data.test = gwas.summary.data.test %>% 
  select(CHR,rsID,A1,MAF,beta.ALT,P,Rsq_ave,N_Controls,
         N_Cases) %>% 
  rename(BETA = beta.ALT,
         INFO = Rsq_ave,
         N0 = N_Controls,
         N1 = N_Cases)
write.table(gwas.summary.data.test,
            file = "/data/zhangh24/multi_ethnic/result/breast_cancer/best_eur_gwas_summary_stat",
            row.names = F,
            col.names = T,
            quote=F)
#idx <- which(gwas.summary.data.test$POS%in%freq.infor$POS==F)
