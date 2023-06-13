Args = commandArgs(trailingOnly = T)
i = as.numeric(Args[[1]])
eth <- c("EUR","AFR","AMR","EAS","SAS")

#system(paste0("cd /dcs04/nilanjan/data/23andme/cleaned; zip -r /dcs04/nilanjan/data/hzhang1/multi_ethnic/data/cleaned/",eth[i],".zip  ",eth[i]))
system(paste0("cd /data/zhangh24/multi_ethnic/data/cleaned; unzip ",eth[i],".zip"))
#system(paste0("cd /data/zhangh24/multi_ethnic/result/cleaned/organize_prs/; zip -r TDLD_EB.zip TDLD_EB"))