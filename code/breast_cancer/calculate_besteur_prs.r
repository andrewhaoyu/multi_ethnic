#best eur prs
setwd("/data/zhangh24/multi_ethnic/data/")
trait = c("overall","erpos","erneg")
library(data.table)
library(tidyverse)

prs.snp <- read.csv("/data/zhangh24/multi_ethnic/data/breast_cancer_330SNPs.csv",header=T)
load("./AABC_data/BC_EUR_overall_aligned.rdata")
prs.snp = prs.snp %>% 
  unite("chr.pos",CHR,Position,sep=":",remove=F)