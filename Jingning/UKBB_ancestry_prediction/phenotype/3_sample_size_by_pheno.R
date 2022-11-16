
rm(list=ls())

library(bigreadr)
library(readr)
library(dplyr)

eth_all <- c("AFR","AMR","EAS","EUR","SAS")
trait_all <- c("bmi","height","HDL","LDL","logTG","nonHDL","TC","TG")

res <- tibble()
for (t in 1:length(trait_all)){
  trait <- trait_all[t]

res_t <- matrix(nrow=5,ncol=2)

for (i in 1:length(eth_all)){
  ethnic <- eth_all[i]

  pheno <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/tuning+validation/",ethnic,"_tuning.txt"))
  pheno = pheno[,1:2]
  covar <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/",ethnic,"_tuning.txt"))
  pheno <- left_join(pheno, covar)
  colnames(pheno) = c('id','y','sex','age',paste0('pc',1:10))
  pheno = pheno[complete.cases(pheno$y),]

  res_t[i,1] <- nrow(pheno)


  pheno <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/tuning+validation/",ethnic,"_validation.txt"))
  pheno = pheno[,1:2]
  covar <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/",ethnic,"_validation.txt"))
  pheno <- left_join(pheno, covar)
  colnames(pheno) = c('id','y','sex','age',paste0('pc',1:10))
  pheno = pheno[complete.cases(pheno$y),]

  res_t[i,2] <- nrow(pheno)
}


  res_t <- data.frame(ethnic=eth_all,trait=trait, res_t)
  colnames(res_t) <- c("ethnic","trait","tuning","validation")

  res <- rbind(res, res_t)
}
write_tsv(res, "/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/complete_case_sample_size.txt")

