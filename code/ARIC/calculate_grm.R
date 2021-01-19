args = commandArgs(trailingOnly = T)
#i is the ethnic
#j is the number of chromsome
i = as.numeric(args[[1]])
j = as.numeric(args[[2]])

eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")

setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic")
temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[1],"/",eth[i],"/")
data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
res <- system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/gcta/gcta64 --bfile ",data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",j," --maf 0.01 --make-grm --out ",temp.dir,"chr.qc",j))
