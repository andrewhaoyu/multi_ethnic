eth_vec <- c("EUR","AFR","AMR")
trait_vec <-c("height","bmi")
library(data.table)
result_list = list()
temp = 1
for(l in 1:2){
  for(i in 2:3){
    eth = eth_vec[i]
    trait = trait_vec[l]
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/AOU/PRSCSX/",eth_vec[i],"/",trait_vec[l],"/")
    load(paste0(out.dir, "prscsx.result"))
    
    result_list[[temp]] = r2.result
    temp = temp+ 1
  }
}

final_result = rbindlist(result_list)
save(final_result, file = "/data/zhangh24/multi_ethnic/result/AOU/PRSCSX/prscsx.rdata")
