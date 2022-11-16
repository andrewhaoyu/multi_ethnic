
rm(list=ls())

library(bigreadr)
library(readr)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

phenotype <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/",trait,"_from_ukbb.txt"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/tuning+validation/"))

for (ethnic in c("AFR","AMR","EAS","EUR","SAS")){

  id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/all_data/",ethnic,"/chr1.fam"))
  write_tsv(phenotype[match(id$V1, phenotype$`f.eid`),], paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/tuning+validation/",ethnic,"_all_data.txt"))

  #tuning_id <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/",ethnic,".tuning_id"))[[1]]
  #write_tsv(phenotype[phenotype$`f.eid` %in% tuning_id, ], paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/tuning+validation/",ethnic,"_tuning.txt"))
  #tuning_id <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/",ethnic,".validation_id"))[[1]]
  #write_tsv(phenotype[phenotype$`f.eid` %in% tuning_id, ], paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/tuning+validation/",ethnic,"_validation.txt"))

}
