library(dplyr)
library(data.table)
setwd("/dcs04/nilanjan/data/rzhao/AllofUs/GWAS/")
eth = c("eur","afr","amr")

for(i in 1:3){
  data.list = list()
  temp = 1
  for(j in 1:22){
    data = fread(paste0("./plink_height/plink_height_",eth[i],"/GWAS_height_",eth[i],"_chr",j,".Height.glm.linear"))  
    data.list[[j]] = data
  }
  data.com = rbindlist(data.list)
  write.table(data.com,
              file = paste0("/dcs04/nilanjan/data/hzhang1/multi_ethnic/data/AOU/height/sum_stat_",eth[i],".txt"),
              quote = F,
              row.names = F)
  }



