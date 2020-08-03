#run GWAS using plink2
#use plink2 to run
args  = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
i = as.numeric(args[[1]])
j = as.numeric(args[[2]])
l = as.numeric(args[[3]])
#load the phenotpypes data and use plink to run
eth <- c("EUR","AFR","AMR","EAS","SAS")
library(data.table)
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"

system(paste0("/data/zhangh24/software/plink2 --bfile ",cur.dir,eth[i],"/chr",j,".tag --out ",cur.dir,eth[i],"/summary_chr",j,"_rho_",l,".out --linear --all-pheno --allow-no-sex --pheno ",cur.dir,eth[i],"/pheno_plink_rho",l))


# for(i in 1:5){
#   fam <- data.frame(fread(paste0(cur.dir,eth[i],"/all_chr.tag.fam")))
#     for(j in 1:22){
#     write.table(fam, file = paste0(cur.dir,eth[i],"/chr",j,".tag.fam"),row.names = F,col.names = F,quote=F)
#   }
# }
# 
# 
