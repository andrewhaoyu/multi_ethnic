args = commandArgs(trailingOnly = T)
i = as.numric(args[[1]])
j = as.numeric(args[[2]])

setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic")
temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
res <- system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/gcta/gcta64 --bfile ",data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",j," --maf 0.01 --make-grm --out ",temp.dir,"chr.qc",j))
