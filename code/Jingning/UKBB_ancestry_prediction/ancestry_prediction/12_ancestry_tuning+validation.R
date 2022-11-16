
rm(list=ls())

library(bigreadr)
library(readr)

exclude <- readLines("/dcl01/arking/data/UK_Biobank/static/Phenotype/Samples.to.Exclude/w17731_20220222.csv")

set.seed(20110715)

for(ethnic in c("AFR","AMR","EAS","SAS")){
  tmp <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/",ethnic,".id"), col_name=F, col_type="cc")
  tmp <- tmp[!(tmp$X2 %in% exclude),]
  #tuning_ind <- sample(1:nrow(tmp), floor(nrow(tmp)/2))
  #validatation_ind <- (1:nrow(tmp))[-tuning_ind]
  #write_tsv(tmp[tuning_ind,], paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/",ethnic,".tuning_id"))
  #write_tsv(tmp[validatation_ind,], paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/",ethnic,".validation_id"))
  write_tsv(tmp, paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/",ethnic,".all_id"))
}


## EUR (unrelated_whites.id from /dcl01/arking/data/UK_Biobank/active/het_homo_mortality/ukb.unrelated.whites.txt)
EUR <- read_tsv("/dcl01/arking/data/UK_Biobank/active/het_homo_mortality/ukb.unrelated.whites.txt")
tmp <- format(EUR$FID, trim = T, scientific = F); tmp <- tmp[!(tmp %in% exclude)]
tmp <- data.frame(FID=tmp, IID=tmp)
#tuning_ind <- sample(1:nrow(tmp), 10000)
#validatation_ind <- sample((1:nrow(tmp))[-tuning_ind], 10000)
#write_tsv(tmp[tuning_ind,], "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/EUR.tuning_id")
#write_tsv(tmp[validatation_ind,], "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/EUR.validation_id")
tuning_ind <- read_tsv("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/EUR.tuning_id")$IID
validatation_ind <- read_tsv("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/EUR.validation_id")$IID
write_tsv(tmp[tmp$IID %in% c(tuning_ind,validatation_ind),], paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/EUR.all_id"))

