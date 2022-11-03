eth = c("EUR","AFR","AMR","EAS","SAS")
trait <- c("HDL","LDL",
           "logTG",
           "TC")
trait_name = c("High-density lipoprotein cholesterol","Low-density lipoprotein cholesterol",
               "Log triglycerides",
               "Total triglycerides")
eth_name = c("European","African American",
             "Latino","East Asian",
             "South Asian")
library(data.table)
total = 4
h2_vec = rep(0,total)
se_vec = rep(0,total)
trait_vec = rep("c",total)
eth_vec  = rep("c",total)
temp=1
for(i in 1:1){
  for(l in 1:4){
    result.folder = paste0("/data/zhangh24/multi_ethnic/result/GLGC/herit/",eth[i],"/",trait[l])
    setwd(result.folder)
    file = readLines("result.log")
    idx = grep("Total Observed scale h2:",file)
    tempstring = gsub("Total Observed scale h2:","",file[idx])
    split.string = strsplit(tempstring,split="\\(")[[1]]
    h2 = as.numeric(split.string[1])
    se = as.numeric(gsub("\\)","",split.string[2]))
    h2_vec[temp] = h2
    se_vec[temp] = se
    trait_vec[temp] = trait_name[l]
    eth_vec[temp] = eth_name[i]
    temp = temp  + 1
  }
}
sigma2toauc = function(x){
  ifelse(x==0,0.50,round(pnorm(sqrt(0.5*x)),3))
}
result.long.h2 = data.frame(eth_vec,trait_vec,   h2_new = paste0(round(h2_vec,3)," (",round(se_vec,3),")"))

library(tidyr)
result.wide = spread(result.long.h2,eth_vec,h2_new)
result.wide
result.long.auc = data.frame(eth_vec,trait_vec,    auc = sigma2toauc(h2_vec))
result.wide.auc = spread(result.long.auc,eth_vec,auc)
result.wide.auc
write.csv(result.wide,"/data/zhangh24/multi_ethnic/result/GLGC/herit/herit.table.csv")
write.csv(result.wide.auc,"/data/zhangh24/multi_ethnic/result/cleaned/herit/herit.table.auc.csv")
