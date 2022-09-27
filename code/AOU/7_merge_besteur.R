eth_vec <- c("EUR","AFR","AMR")
trait_vec <- c("height","bmi")
library(data.table)
result = list()
temp = 1
for(l in 1:2){
  for(i in 2:3){
    eth = eth_vec[i]
    trait = trait_vec[l]
    out.dir = out.dir = paste0("/data/zhangh24/multi_ethnic/result/AOU/BestEUR/",eth,"/",trait,"/")
    load(paste0(out.dir, "BestEUR.result"))
    result[[temp]] = r2.result
    temp = temp+ 1
  }
}

final_result = rbindlist(result)
save(final_result, file = "/data/zhangh24/multi_ethnic/result/AOU/BestEUR/best_eur.rdata")
