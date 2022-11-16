
rm(list=ls())

library(bigreadr)
library(readr)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/nonHDL"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/nonHDL/tuning+validation/"))

for(ancestry in c("AFR","AMR","EAS","EUR","SAS")){
  TC <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/TC/tuning+validation/",ancestry,"_all_data.txt"), col_types=cols())
  HDL <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/HDL/tuning+validation/",ancestry,"_all_data.txt"), col_types=cols())
  if(identical(TC[[1]], HDL[[1]])){
    print(paste0(ancestry, ": sample identical!"))
    nonHDL <- data.frame("f.eid" = TC$f.eid, "f.30690.0.0minusf.30760.0.0" = TC[[2]] - HDL[[2]])
    write_tsv(nonHDL, paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/nonHDL/tuning+validation/",ancestry,"_all_data.txt"))
  }else{
    print(paste0(ancestry, ": sample NOT identical! SKIPPED!"))
  }

}

