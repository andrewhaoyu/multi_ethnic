setwd("/data/zhangh24/multi_ethnic/data/")
load("./AABC_data/BC_EUR_overall_aligned.rdata")
library(tidyverse)
library(data.table)
ukb.sum <- fread("ukb_gera_breast_meta.txt",header=T)
#add a little bit noise if OR is 1
#otherwise, se can't be calculated based on P value

ukb.sum = ukb.sum %>% 
  mutate(OR_update = ifelse(OR==1,1.0001,OR))

ukb.sum.update = ukb.sum %>% 
  unite("chr.pos",CHR,BP,sep = ":",remove=F) %>% 
  rename(UKB_effect_allele = A1,
         UKB_other_allele = A2,
         POS=BP,
         rsid = ID) %>% 
  mutate(beta_ukb=log(OR_update),
         absz = qnorm(P/2),
         se_ukb = abs(beta_ukb/absz))
ukb.sum.update = ukb.sum.update %>% 
  select(chr.pos,UKB_effect_allele,
         UKB_other_allele,
         beta_ukb,se_ukb,rsid)

# idx <- which(sum.data$ID=="1:100880328:A:T")
# sum.data[idx,]
# jdx <- which(ukb.sum.update$chr.pos=="4:84370124")
# ukb.sum.update[jdx,]

#match the ukb SNPs with BCAC SNPs
#align alleles
#keep the BCAC only SNPs
#keep two multi-allelic EUR SNPs 
head(ukb.sum.update)
idx <- which(ukb.sum.update$chr.pos%in%
               c("4:84370124","7:139943702"))
ukb.sum.update$beta_ukb[idx] = NA

bc_aligned = left_join(sum.data,
                       ukb.sum.update,by="chr.pos") %>% 
  filter(is.na(beta_ukb)|(((Eff_allele==UKB_effect_allele)&(Ref_allele==UKB_other_allele))|
  ((Eff_allele==UKB_other_allele)&(Ref_allele==UKB_effect_allele))))



bc_aligned = bc_aligned %>% 
  mutate(beta_ukb_update = ifelse(UKB_effect_allele==Eff_allele,beta_ukb,-beta_ukb))

Meta = function(coef_vec,se_vec){
  var_vec = se_vec^2
  meta_var = (sum(1/var_vec))^-1
  meta_coef = meta_var*sum(coef_vec/var_vec)
  meta_se = sqrt(meta_var)
  return(c(meta_coef,meta_se))
}
N = nrow(bc_aligned)
beta_meta = rep(0,N)
se_meta = rep(0,N)
Beta_eur_update = bc_aligned$Beta_eur_update
se_eur = bc_aligned$Se_eur

ukb_beta = bc_aligned$beta_ukb_update
ukb_se = bc_aligned$se_ukb

for(k in 1:N){
  if(k%%1000==0){print(k)}
  coef_vec = c(Beta_eur_update[k],ukb_beta[k])
  se_vec = c(se_eur[k],ukb_se[k])
  if(is.na(coef_vec[2])){
    beta_meta[k] = coef_vec[1]
    se_meta[k] = se_vec[1]
  }else{
    meta_result = Meta(coef_vec,se_vec)
    beta_meta[k] = meta_result[1]
    se_meta[k] = meta_result[2]
  }
  
}
bc_aligned$beta_meta = beta_meta
bc_aligned$se_meta = se_meta
bc_aligned$p_meta = 2*pnorm(-abs(beta_meta/se_meta),lower.tail = T)
bc_aligned_update = bc_aligned %>% 
  select(ID,chr.pos,CHR,POS,Eff_allele,Ref_allele,
         beta_meta,se_meta,p_meta,
         R2_eur,
         MAF_eur,KG.ID,rsid) %>% 
  rename(Beta_eur_update = beta_meta,
         Se_eur = se_meta,
         p_eur = p_meta)
sum.data = bc_aligned_update
#idx <- which(bc_aligned$chr.pos=="4:84370124")
save(sum.data, file = paste0("./AABC_data/BC_EUR_overall_ukb_meta_aligned.rdata"))
