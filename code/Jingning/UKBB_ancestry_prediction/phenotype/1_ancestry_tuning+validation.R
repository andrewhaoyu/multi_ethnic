
rm(list=ls())

library(bigreadr)
library(readr)

ancestry <- fread2("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/ethnic_background_from_ukbb.txt")
exclude <- readLines("/dcl01/arking/data/UK_Biobank/static/Phenotype/Samples.to.Exclude/w17731_20220222.csv")

write_tsv(data.frame(FID=exclude, IID=exclude), "/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/tuning+validation/exclude.id")


set.seed(20110711)

## SAS (Indian, Pakistani, Bangladeshi)
SAS <- ancestry[ancestry$`f.21000.0.0` %in% c(3001,3002,3003),]
tmp <- format(SAS$f.eid, trim = T, scientific = F); tmp <- tmp[!(tmp %in% exclude)]
tmp <- data.frame(FID=tmp, IID=tmp)
tuning_ind <- sample(1:nrow(tmp), floor(nrow(tmp)/2))
validatation_ind <- (1:nrow(tmp))[-tuning_ind]
write_tsv(tmp[tuning_ind,], "/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/tuning+validation/SAS.tuning_id")
write_tsv(tmp[validatation_ind,], "/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/tuning+validation/SAS.validation_id")


## EAS (Chinese)
EAS <- ancestry[ancestry$`f.21000.0.0` %in% c(5),]
tmp <- format(EAS$f.eid, trim = T, scientific = F); tmp <- tmp[!(tmp %in% exclude)]
tmp <- data.frame(FID=tmp, IID=tmp)
tuning_ind <- sample(1:nrow(tmp), floor(nrow(tmp)/2))
validatation_ind <- (1:nrow(tmp))[-tuning_ind]
write_tsv(tmp[tuning_ind,], "/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/tuning+validation/EAS.tuning_id")
write_tsv(tmp[validatation_ind,], "/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/tuning+validation/EAS.validation_id")


## AFR (Caribbean, African)
AFR <- ancestry[ancestry$`f.21000.0.0` %in% c(4001,4002),]
tmp <- format(AFR$f.eid, trim = T, scientific = F); tmp <- tmp[!(tmp %in% exclude)]
tmp <- data.frame(FID=tmp, IID=tmp)
tuning_ind <- sample(1:nrow(tmp), floor(nrow(tmp)/2))
validatation_ind <- (1:nrow(tmp))[-tuning_ind]
write_tsv(tmp[tuning_ind,], "/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/tuning+validation/AFR.tuning_id")
write_tsv(tmp[validatation_ind,], "/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/tuning+validation/AFR.validation_id")

## EUR (unrelated_whites.id from /dcl01/arking/data/UK_Biobank/active/het_homo_mortality/ukb.unrelated.whites.txt)
EUR <- read_tsv("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/unrelated_whites.id")
tmp <- format(EUR$FID, trim = T, scientific = F); tmp <- tmp[!(tmp %in% exclude)]
tmp <- data.frame(FID=tmp, IID=tmp)
tuning_ind <- sample(1:nrow(tmp), 10000)
validatation_ind <- sample((1:nrow(tmp))[-tuning_ind], 10000)
write_tsv(tmp[tuning_ind,], "/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/tuning+validation/EUR.tuning_id")
write_tsv(tmp[validatation_ind,], "/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/tuning+validation/EUR.validation_id")


