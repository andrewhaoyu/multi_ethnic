eth_vec <- c("EUR","AFR","AMR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
library(data.table)
result_list = list()
temp = 1
for(l in 1:4){
  for(i in 1:5){
    eth = eth_vec[i]
    trait = trait_vec[l]
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/",eth,"/",trait,"/")
    load(paste0(out.dir, "CTSLEB_all.result"))
    
    
    result_list[[temp]] = result
    temp = temp+ 1
  }
}

final_result = rbindlist(result_list)
save(final_result, file = "/data/zhangh24/multi_ethnic/result/GLGC/CTSLEB/ct_sleb_all.rdata")



#test results

eth_vec <- c("EUR","AFR","AMR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
library(data.table)
result_list = list()
temp = 1
for(l in 1:4){
  for(i in c(2,4,5)){
    eth = eth_vec[i]
    trait = trait_vec[l]
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/",eth,"/",trait,"/")
    load(paste0(out.dir, "CTSLEB_all.result"))
    
    
    result_list[[temp]] = result
    temp = temp+ 1
  }
}

final_result = rbindlist(result_list)

eth_vec <- c("EUR","AFR","AMR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
library(data.table)
result_list = list()
temp = 1
for(l in 1:4){
  for(i in c(2,4,5)){
    eth = eth_vec[i]
    trait = trait_vec[l]
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/",eth,"/",trait,"/")
    load(paste0(out.dir, "CTSLEB_all_pgs.result"))
    
    
    result_list[[temp]] = data.frame(eth = eth,trait = trait,
      r2 = r2_ctsleb)
    temp = temp+ 1
  }
}

final_result2 = rbindlist(result_list)
head(final_result)
head(final_result2)
final_result$r2/final_result2$r2
