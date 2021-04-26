#best eur prs
setwd("/data/zhangh24/multi_ethnic/data/")
trait = c("overall","erpos","erneg")
library(data.table)
library(tidyverse)

prs.snp <- read.csv("/data/zhangh24/multi_ethnic/data/breast_cancer_330SNPs.csv",header=T)
prs.snp = prs.snp %>% 
  select(variant,chr.pos,Reference.allele,Effect.allele)
load(paste0("./AABC_data/BC_EUR_",trait[l],"_aligned.rdata"))
prs.snp = prs.snp %>% 
  unite("chr.pos",CHR,Position,sep=":",remove=F)
prs.snp.eur = left_join(prs.snp,sum.data,by="chr.pos") %>% 
  filter(((Eff_allele==Effect.allele)&(Ref_allele==Reference.allele))|
           (Eff_allele==Reference.allele)&(Ref_allele==Effect.allele))

# jdx <- which(sum.data$chr.pos=="7:91459189")
# sum.data[jdx,]
idx <- which(prs.snp$variant%in%prs.snp.eur$variant==F)
prs.snp[idx[1:10],]
