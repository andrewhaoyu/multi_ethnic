
rm(list=ls())

library(bigreadr)
library(readr)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/logTG"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/logTG/tuning+validation/"))

for(ancestry in c("AFR","AMR","EAS","EUR","SAS")){

  TG <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/TG/tuning+validation/",ancestry,"_all_data.txt"), col_types=cols())
  logTG <- data.frame("f.eid" = TG$f.eid, "log.f.30870.0.0" = log(TG$f.30870.0.0))
  write_tsv(logTG, paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/logTG/tuning+validation/",ancestry,"_all_data.txt"))

}

