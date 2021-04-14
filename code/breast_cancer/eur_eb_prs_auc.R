setwd("/data/zhangh24/multi_ethnic/data/AABC_data")
load("BC_AFR_overall_train_KGID.rdata")
sum.data.train = sum.data.update
library(data.table)
library(dplyr)
prs.snp <- read.csv("/data/zhangh24/multi_ethnic/data/breast_cancer_330SNPs.csv",header=T)
#load eur data
sum.eur <- fread("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/prepare_summary_level_statistics/result/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt",header=T)
colnames(sum.eur)[1] <- "Var_name"

#idx <- which(sum.eur$Position.Onco==29121087)


prs.snp.eur <- left_join(prs.snp,sum.eur,by="Var_name") 

prs.eur = prs.snp.eur %>%  
  select(Var_name,variant,CHR,Position,Effect.allele,
         EAF3,Beta.meta,sdE.meta,p.meta) %>% 
  rename(POS = Position,
         Effect_allele_eur=Effect.allele,
         EAF.eur = EAF3,
         Beta.eur = Beta.meta,
         se.eur = sdE.meta,
         p.eur = p.meta,
         rsid=variant)

idx <- which(prs.eur$Var_name=="")


var_name = data.frame(ID = gsub("_",":",prs.snp[,1]),rsID = prs.snp[,2])
# idx <- which(var_name%in%sum.data[,1]==F)
# 
# jdx <- which(sum.data$CHR==17&
#                sum.data$POS==7571752)
library(dplyr)
prs.tar.all = inner_join(var_name,sum.data.train,
                          by="ID")

prs.tar = prs.tar.all %>% 
  rename(Effect_allele_tar = ALT,
         Beta.tar = beta.ALT,
         se.tar = SE,
         p.tar = P,
         rsid = rsID) %>% 
  mutate(MAF = ifelse(AF_ALT<=0.5,AF_ALT,1-AF_ALT))%>% 
  select(rsid,Effect_allele_tar,
         EAF.tar = MAF,
         Beta.tar,se.tar,p.tar)
  
#combine the two ethnic groups
prs.all <- left_join(prs.tar,
                     prs.eur,
                     by="rsid") %>% 
  mutate(Beta.eur = ifelse((Effect_allele_eur==EAF.tar)|
                             is.na(EAF.tar),
                           Beta.eur,-Beta.eur))
# idx <- which(prs.all$Effect_allele_eur!=
#                prs.all$Effect_allele_tar)
all.equal(prs.all$Effect_allele_eur,prs.all$Effect_allele_tar)





beta_tar <- summary.com.prior$beta_tar
sd_tar <- summary.com.prior$sd_tar
beta_eur <- summary.com.prior$beta_eur
sd_eur <- summary.com.prior$sd_eur

EBprior = EstimatePrior(beta_tar,sd_tar,
                        beta_eur,sd_eur)







#prepare data for summary AUC
prs.model.file = best.eur.snp %>% 
  mutate(Effect_Allele = ALT,
         BETA = beta.ALT) %>% 
  select(KG.ID,Effect_Allele,BETA)

write.table(prs.model.file,
            file = "/data/zhangh24/multi_ethnic/result/breast_cancer/best_eur_prsmodelfile",
            row.names = F,
            col.names = F,
            quote=F,
            sep = "\t")
load("/data/zhangh24/multi_ethnic/data/AABC_data/BC_AFR_overall_valid_KGMAF.rdata")
sum.data.vad = sum.data.update
best.eur.snp.select = best.eur.snp %>% 
  select(ID,KG.ID,ALT,MAF) %>% 
  rename(A1 = ALT)

gwas.summary.data.test = left_join(best.eur.snp.select,
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
