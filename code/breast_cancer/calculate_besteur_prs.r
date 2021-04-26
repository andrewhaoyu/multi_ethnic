#best eur prs
setwd("/data/zhangh24/multi_ethnic/data/")
#load GHBS_id
bim <- fread("/data/zhangh24/multi_ethnic/data/GBHS_plink/all_chr_temp.bim",header=F) 
bim = bim %>% 
  unite("chr.pos",V1,V4,sep=":",remove=F)
bim = bim %>% 
  rename(GBHS_id = V2)


trait = c("overall","erpos","erneg")
library(data.table)
library(tidyverse)

prs.snp <- read.csv("/data/zhangh24/multi_ethnic/data/breast_cancer_330SNPs.csv",header=T)
prs.snp = prs.snp %>% 
  unite("chr.pos",CHR,Position,sep=":",remove=F)

prs.snp = prs.snp %>% 
  select(variant,chr.pos,Reference.allele,Effect.allele)
load(paste0("./AABC_data/BC_EUR_",trait[l],"_aligned.rdata"))
prs.snp.eur = left_join(prs.snp,sum.data,by="chr.pos") %>% 
  filter(((Eff_allele==Effect.allele)&(Ref_allele==Reference.allele))|
           (Eff_allele==Reference.allele)&(Ref_allele==Effect.allele)) %>% 
  select(c(colnames(sum.data),"chr.pos","variant"))

jdx <- which(bim$chr.pos=="1:114445880")



prs.snp.eur.update = left_join(prs.snp.eur,bim,by="chr.pos") %>% 
  filter(((Eff_allele==V5)&(Ref_allele==V6))|
           (Eff_allele==V6)&(Ref_allele==V5))

idx <- which(prs.snp.eur$ID%in%prs.snp.eur.update$ID==F)
