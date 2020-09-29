#
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
ethnic <- c("EUR","AFR","AMR","EAS","SAS")
system(paste0("cd /data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i1]))
system(paste0("summary_out_rep"))
system(paste0("cp summary_out_rho_*_rep_* summary_out_rep/"))
system(paste0("tar -zcvf summary_out_rep.tar.gz summary_out_rep"))