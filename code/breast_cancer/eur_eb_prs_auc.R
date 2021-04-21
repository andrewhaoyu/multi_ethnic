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





#chr.pos = paste0(prs.eur$CHR,":",prs.eur$POS)

#idx <- which(temp%in%chr.pos)

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
  mutate(Beta.eur = ifelse((Effect_allele_eur==Effect_allele_tar)|
                             is.na(Effect_allele_tar),
                           Beta.eur,-Beta.eur))
# idx <- which(prs.all$Effect_allele_eur!=
#                prs.all$Effect_allele_tar)
all.equal(prs.all$Effect_allele_eur,prs.all$Effect_allele_tar)


source("/data/zhangh24/multi_ethnic/code/stratch/EB_function.R")


beta_tar <- prs.all$Beta.tar
sd_tar <- prs.all$se.tar
beta_eur <- prs.all$Beta.eur
sd_eur <- prs.all$se.eur

EBprior = EstimatePrior(beta_tar,sd_tar,
                        beta_eur,sd_eur)

post_beta_mat = EBpost(beta_tar,sd_tar,beta_eur,sd_eur,EBprior)
post_beta_tar = post_beta_mat[,1,drop=F]

prs.tar.all$BETA = post_beta_tar

#prepare data for summary AUC
prs.model.file = prs.tar.all %>% 
  mutate(Effect_Allele = ALT) %>% 
  select(KG.ID,Effect_Allele,BETA)

write.table(prs.model.file,
            file = "/data/zhangh24/multi_ethnic/result/breast_cancer/best_eur_eb_prsmodelfile",
            row.names = F,
            col.names = F,
            quote=F,
            sep = "\t")
load("/data/zhangh24/multi_ethnic/data/AABC_data/BC_AFR_overall_valid_KGID.rdata")
sum.data.vad = sum.data.update
best.eur.snp.select = prs.tar.all %>% 
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
source("/data/zhangh24/multi_ethnic/code/breast_cancer/summaryAUC.R")
res = auc(prs.model.file = "/data/zhangh24/multi_ethnic/result/breast_cancer/best_eur_eb_prsmodelfile", 
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

