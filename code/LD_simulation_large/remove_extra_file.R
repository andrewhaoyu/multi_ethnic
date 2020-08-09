#remove extra file to clean space
#i represent ehtnic group
#j represent chr
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
j = as.numeric(args[[2]])
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
system(paste0("rm ",cur.dir,eth[i],"/chr",j,"_*",".tag.*"))
system(paste0("rm ",cur.dir,eth[i],"/summary_chr",j,"*"))
