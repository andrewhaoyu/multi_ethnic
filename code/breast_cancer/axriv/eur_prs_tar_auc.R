setwd("/data/zhangh24/multi_ethnic/data/AABC_data")
load("BC_AFR_overall_train_KGID.rdata")
sum.data.train = sum.data.update
prs.snp <- read.csv("/data/zhangh24/multi_ethnic/data/breast_cancer_330SNPs.csv",header=T)

var_name = data.frame(ID = gsub("_",":",prs.snp[,1]),rsID = prs.snp[,2])
# idx <- which(var_name%in%sum.data[,1]==F)
# 
# jdx <- which(sum.data$CHR==17&
#                sum.data$POS==7571752)
library(dplyr)
best.eur.snp = inner_join(var_name,sum.data.train,
                         by="ID")








# best.eur.snp = best.eur.snp %>%
#   mutate(chr_pos = paste0(CHR,":",POS))
# load("/data/zhangh24/KG.plink/AFR/chr_all_frq.rdata")
# freq.infor = freq.infor %>%
#   mutate(chr_pos = paste0(CHR,":",POS))
# freq.infor.sub = freq.infor %>%
#   rename(KG.ID = SNP) %>%
#   select(KG.ID,
#          chr_pos,MAF)
# 
# best.eur.snp = inner_join(best.eur.snp,
#                           freq.infor.sub,
#                           by="chr_pos")


#prepare data for summary AUC
prs.model.file = best.eur.snp %>% 
  mutate(Effect_Allele = ALT,
         BETA = beta.ALT) %>% 
  select(KG.ID,Effect_Allele,BETA)

write.table(prs.model.file,
            file = "/data/zhangh24/multi_ethnic/result/breast_cancer/best_eur_tarcoef_prsmodelfile",
            row.names = F,
            col.names = F,
            quote=F,
            sep = "\t")
load("BC_AFR_overall_valid_KGID.rdata")
sum.data.vad = sum.data.update
best.eur.snp.select = best.eur.snp %>% 
  mutate(MAF = ifelse(AF_ALT<=0.5,AF_ALT,1-AF_ALT)) %>% 
  select(ID,ALT,MAF) %>% 
  rename(A1 = ALT)

gwas.summary.data.test = inner_join(best.eur.snp.select,
                                   sum.data.vad,by="ID")
all.equal(gwas.summary.data.test$A1,
          gwas.summary.data.test$ALT)
library(data.table)
#gwas.summary.data.test = left_join(gwas.summary.data.test,freq.infor,by="POS")


gwas.summary.data.test = gwas.summary.data.test %>% 
  select(CHR,KG.ID,A1,MAF,beta.ALT,P,Rsq_ave) %>% 
  rename(BETA = beta.ALT,
         INFO = Rsq_ave,
         SNP = KG.ID)
write.table(gwas.summary.data.test,
            file = "/data/zhangh24/multi_ethnic/result/breast_cancer/best_eur_gwas_summary_stat",
            row.names = F,
            col.names = T,
            quote=F)
#idx <- which(gwas.summary.data.test$POS%in%freq.infor$POS==F)
res = auc(prs.model.file = "/data/zhangh24/multi_ethnic/result/breast_cancer/best_eur_tarcoef_prsmodelfile", 
          gwas.summary.stats.file = "/data/zhangh24/multi_ethnic/result/breast_cancer/best_eur_gwas_summary_stat",
          N0 = 3119,
          N1 = 2702,
          soFile = '/data/zhangh24/multi_ethnic/code/breast_cancer/getAdjCorrelation.so',
          flag.correlation.adj.imputated.data = FALSE,
          pos_thr = 5e8,
          KG.plink.pre = '/data/zhangh24/KG.plink/AFR/chr_all')
cat("\n#######################################\n\n")

cat("Predicted AUC:\t", res[1], "\n", sep = "")
cat("Predicted AUC's variance:\t", res[2], "\n", sep = "")

cat("\nHave a nice Day!\n")

cat("\n#######################################\n")