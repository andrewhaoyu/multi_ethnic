setwd("/data/zhangh24/multi_ethnic/data/")
trait = c("overall","erpos","erneg")
library(data.table)
library(tidyverse)
#process overall bc
#load AFR coefficients
l =1 
load(paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS.rdata"))
sum.data = sum.data %>% 
  unite("chr.pos",CHR,POS,sep=":",remove=F)
#idx <- which(sum.data$chr.pos=="1:121287994")
#load eur coefficients
sum.eur <- fread("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/prepare_summary_level_statistics/result/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt",header=T)
sum.eur.update = sum.eur %>% 
  select(var_name,
         Beta.meta,
         sdE.meta,p.meta,
         chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,
         EAFcontrols.Onco,r2.Onco) %>% 
  rename(CHR = chr.Onco,
         POS = Position.Onco,
         Eff_allele_eur = Effect.Onco,
         Ref_allele_eur = Baseline.Onco,
         EAF_eur = EAFcontrols.Onco,
         R2_eur = r2.Onco,
         Beta_eur = Beta.meta,
         Se_eur = sdE.meta,
         p_eur = p.meta)
sum.eur.update = sum.eur.update %>% 
  unite("chr.pos",CHR,POS,sep=":",remove=T)

sum.eur.match = inner_join(sum.data,
                           sum.eur.update,
                           by="chr.pos") %>% 
  filter(((Eff_allele==Eff_allele_eur)&(Ref_allele==Ref_allele_eur))|
           ((Eff_allele==Ref_allele_eur)&(Ref_allele==Eff_allele_eur)))

sum.eur.match.update = sum.eur.match %>% 
  mutate(Beta_eur_update = ifelse(Eff_allele_eur==Eff_allele,Beta_eur,-Beta_eur),
         MAF_eur = ifelse(EAF_eur<=0.5,EAF_eur,1-EAF_eur)) %>% 
  select(ID,chr.pos,CHR,POS,Eff_allele,Ref_allele,Beta_eur_update,
         Se_eur,p_eur,R2_eur,MAF_eur,KG.ID) 
sum.data = sum.eur.match.update
save(sum.data,file ="./AABC_data/BC_EUR_overall_aligned.rdata")






#process ER+ bc
#load AFR coefficients
l =2
load(paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS.rdata"))
sum.data = sum.data %>% 
  unite("chr.pos",CHR,POS,sep=":",remove=F)
#load eur coefficients
sum.eur <- fread("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/prepare_summary_level_statistics/result/oncoarray_bcac_public_release_oct17.txt",header=T)
#colnames(sum.eur)[1] <- "Var_name"
sum.eur.update = sum.eur %>% 
  select(var_name,
         bcac_onco_icogs_gwas_erpos_beta,
         bcac_onco_icogs_gwas_erpos_se,bcac_onco_icogs_gwas_erpos_P1df,
         chr,position_b37,a1,a0,
         bcac_onco_icogs_gwas_erpos_eaf_controls,bcac_onco2_erpos_r2) %>% 
  rename(CHR = chr,
         POS = position_b37,
         Eff_allele_eur = a1,
         Ref_allele_eur = a0,
         EAF_eur = bcac_onco_icogs_gwas_erpos_eaf_controls,
         R2_eur = bcac_onco2_erpos_r2,
         Beta_eur = bcac_onco_icogs_gwas_erpos_beta,
         Se_eur = bcac_onco_icogs_gwas_erpos_se,
         p_eur = bcac_onco_icogs_gwas_erpos_P1df)
sum.eur.update = sum.eur.update %>% 
  unite("chr.pos",CHR,POS,sep=":",remove=T)

sum.eur.match = inner_join(sum.data,
                           sum.eur.update,
                           by="chr.pos") %>% 
  filter(((Eff_allele==Eff_allele_eur)&(Ref_allele==Ref_allele_eur))|
           ((Eff_allele==Ref_allele_eur)&(Ref_allele==Eff_allele_eur)))


# Beta_eur = as.numeric(sum.eur.match$Beta_eur)
# idx <- which(is.na(Beta_eur))
# sum.eur.match[idx[1:10],]
# 
# jdx <- which(sum.eur$var_name=="1_855949_G_A")
# sum.eur[jdx,]

sum.eur.match.update = sum.eur.match %>% 
  mutate(Beta_eur =as.numeric(Beta_eur)) %>% 
  filter(!is.na(Beta_eur)) %>% 
  mutate(Se_eur = as.numeric(Se_eur),
         p_eur = as.numeric(p_eur),
         EAF_eur = as.numeric(EAF_eur),
         R2_eur = as.numeric(EAF_eur),
         Beta_eur_update = ifelse(Eff_allele_eur==Eff_allele,Beta_eur,-Beta_eur),
         MAF_eur = ifelse(EAF_eur<=0.5,EAF_eur,1-EAF_eur)) %>% 
  
  select(ID,chr.pos,CHR,POS,Eff_allele,Ref_allele,Beta_eur_update,
         Se_eur,p_eur,R2_eur,MAF_eur,KG.ID) 
sum.data = sum.eur.match.update
#idx <- which(is.na(sum.eur.match.update$Beta_eur_update)
save(sum.data,file ="./AABC_data/BC_EUR_erpos_aligned.rdata")







l =3
load(paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS.rdata"))
sum.data = sum.data %>% 
  unite("chr.pos",CHR,POS,sep=":",remove=F)
#load eur coefficients
#sum.eur <- fread("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/prepare_summary_level_statistics/result/oncoarray_bcac_public_release_oct17.txt",header=T)
#colnames(sum.eur)[1] <- "Var_name"
sum.eur.update = sum.eur %>% 
  select(var_name,
         bcac_onco_icogs_gwas_erneg_beta,
         bcac_onco_icogs_gwas_erneg_se,bcac_onco_icogs_gwas_erneg_P1df,
         chr,position_b37,a1,a0,
         bcac_onco_icogs_gwas_erneg_eaf_controls,bcac_onco2_erneg_r2) %>% 
  rename(CHR = chr,
         POS = position_b37,
         Eff_allele_eur = a1,
         Ref_allele_eur = a0,
         EAF_eur = bcac_onco_icogs_gwas_erneg_eaf_controls,
         R2_eur = bcac_onco2_erneg_r2,
         Beta_eur = bcac_onco_icogs_gwas_erneg_beta,
         Se_eur = bcac_onco_icogs_gwas_erneg_se,
         p_eur = bcac_onco_icogs_gwas_erneg_P1df)
sum.eur.update = sum.eur.update %>% 
  unite("chr.pos",CHR,POS,sep=":",remove=T)

sum.eur.match = inner_join(sum.data,
                           sum.eur.update,
                           by="chr.pos") %>% 
  filter(((Eff_allele==Eff_allele_eur)&(Ref_allele==Ref_allele_eur))|
           ((Eff_allele==Ref_allele_eur)&(Ref_allele==Eff_allele_eur)))


# Beta_eur = as.numeric(sum.eur.match$Beta_eur)
# idx <- which(is.na(Beta_eur))
# sum.eur.match[idx[1:10],]
# 
# jdx <- which(sum.eur$var_name=="1_855949_G_A")
# sum.eur[jdx,]

sum.eur.match.update = sum.eur.match %>% 
  mutate(Beta_eur =as.numeric(Beta_eur)) %>% 
  filter(!is.na(Beta_eur)) %>% 
  mutate(Se_eur = as.numeric(Se_eur),
         p_eur = as.numeric(p_eur),
         EAF_eur = as.numeric(EAF_eur),
         R2_eur = as.numeric(EAF_eur),
         Beta_eur_update = ifelse(Eff_allele_eur==Eff_allele,Beta_eur,-Beta_eur),
         MAF_eur = ifelse(EAF_eur<=0.5,EAF_eur,1-EAF_eur)) %>% 
  
  select(ID,chr.pos,CHR,POS,Eff_allele,Ref_allele,Beta_eur_update,
         Se_eur,p_eur,R2_eur,MAF_eur,KG.ID) 
sum.data = sum.eur.match.update
#idx <- which(is.na(sum.eur.match.update$Beta_eur_update)
save(sum.data,file ="./AABC_data/BC_EUR_erneg_aligned.rdata")
