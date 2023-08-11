eth_vec <- c("EUR","AFR","AMR")
trait_vec <-c("height","bmi")
library(data.table)
result_list = list()
temp = 1
for(l in 1:2){
  for(i in 1:2){
    eth = eth_vec[i]
    trait = trait_vec[l]
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/AOU/clumping_result/PT/",eth,"/",trait,"/")
    load(paste0(out.dir, "CTSLEB_all.result"))
    
    
    result_list[[temp]] = result
    temp = temp+ 1
  }
}

final_result = rbindlist(result_list)
save(final_result, file = "/data/zhangh24/multi_ethnic/result/AOU/CTSLEB/ct_sleb_all.rdata")



eth_vec <- c("EUR","AFR","AMR")
trait_vec <-c("height","bmi")
library(data.table)
result_list = list()
temp = 1
for(l in 1:2){
  for(i in 2:2){
    eth = eth_vec[i]
    trait = trait_vec[l]
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/AOU/clumping_result/PT/",eth,"/",trait,"/")
    load(paste0(out.dir, "CTSLEB_all_pgs.result"))
    
    
    result_list[[temp]] = result
    temp = temp+ 1
  }
}

final_result2 = rbindlist(result_list)
