eth_vec <- c("EUR","AFR","AMR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
library(data.table)
result_list = list()
temp = 1
for(l in 1:2){
  for(i in 2:3){
    eth = eth_vec[i]
    trait = trait_vec[l]
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/polypred/",eth_vec[i],"/",trait_vec[l],"/")
    load(paste0(out.dir, "polypred.result"))
    
    result_list[[temp]] = r2.result
    temp = temp+ 1
  }
}

final_result = rbindlist(result_list)
save(final_result, file = "/data/zhangh24/multi_ethnic/result/GLGC/polypred/polypred.rdata")
