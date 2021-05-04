setwd("/data/zhangh24/multi_ethnic/data/")
load("./AABC_data/BC_EUR_overall_aligned.rdata")
library(tidyverse)
library(data.table)
ukb.sum <- fread("ukb_gera_breast_meta.txt",header=T)

ukb.sum.update = ukb.sum %>% 
  unite("chr.pos",CHR,BP,sep = ":",remove=F) %>% 
  rename(UKB_effect_allele = A1,
         UKB_other_allele = A2,
         POS=BP,
         rsid = ID) %>% 
  mutate(beta_ukb=log(OR),
         absz = qnorm(P/2),
         se_ukb = abs(beta_ukb/absz))
ukb.sum.update = ukb.sum.update %>% 
  select(chr.pos,UKB_effect_allele,
         UKB_other_allele,
         beta_ukb,se_ukb,rsid)
# idx <- which(sum.data$ID=="1:100880328:A:T")
# sum.data[idx,]
# jdx <- which(ukb.sum.update$ID=="rs612683")
# ukb.sum.update[jdx,]

#match the ukb SNPs with BCAC SNPs
#align alleles
#keep the BCAC only SNPs
bc_aligned = left_join(sum.data,
                       ukb.sum.update,by="chr.pos") %>% 
  filter(is.na(beta_ukb)|(((Eff_allele==UKB_effect_allele)&(Ref_allele==UKB_other_allele))|
  ((Eff_allele==UKB_other_allele)&(Ref_allele==UKB_effect_allele))))


idx <- which(sum.data$chr.pos%in%bc_aligned$chr.pos==F)
sum.data[idx[1:10],]
