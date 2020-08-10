#compress phenotypes and summary stat to transfer
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
eth <- c("EUR","AFR","AMR","EAS","SAS")
setwd(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i]))
system(paste0("rm -r pheno_summary_stat"))
system(paste0("mkdir pheno_summary_stat"))
system(paste0("cp summary_out* pheno_summary_stat"))
system(paste0("cp summary_MAF_* pheno_summary_stat"))
system(paste0("cp pheno_plink* pheno_summary_stat"))
system(paste0("cp phenotypes_*.phen pheno_summary_stat"))
system(paste0("tar -zcvf pheno_summary_stat.tar.gz pheno_summary_stat"))
c 