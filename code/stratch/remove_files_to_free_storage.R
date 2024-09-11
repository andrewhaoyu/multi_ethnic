#delete unused files to clean up storage
Args = commandArgs(trailingOnly = T)
i = as.numeric(Args[[1]])

eth = c("EUR","AMR","AFR",
        "EAS","SAS")
system(paste0("cd /data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/;cd ",eth[i],";rm -rf prs"))