
rm(list=ls())

library(bigreadr)

ethnic <- c("EUR","EAS","SAS","AFR","AMR")
M <- length(ethnic)

for (mm in 1:M){
  for (chr in 1:22){
    print(paste0(mm,"-",chr))
    if(chr==1){
      bim <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/all_data/",ethnic[mm],"/chr", chr,".bim"))
    }else{
      bim <- rbind(bim,
                    fread2(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/all_data/",ethnic[mm],"/chr", chr,".bim")))
    }
  }
  fwrite2(bim, paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/all_data/",ethnic[mm],"/allchr_bim"), col.names = F, sep="\t", nThread=1)
}

