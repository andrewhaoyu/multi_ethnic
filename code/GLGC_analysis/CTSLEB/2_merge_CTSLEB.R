eth_vec <- c("EUR","AFR","AMR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
library(data.table)
result_list = list()
temp = 1
for(l in 1:4){
  for(i in 2:5){
    eth = eth_vec[i]
    trait = trait_vec[l]
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/",eth,"/",trait,"/")
    load(paste0(out.dir, "CTSLEB.result"))
    # r2.result = data.frame(eth = eth_vec[i],
    #                        trait = trait_vec[l],
    #                        method = "CT-SLEB",
    #                        r2 = r2_ctsleb
    # )
    
    result_list[[temp]] = result
    temp = temp+ 1
  }
}

final_result = rbindlist(result_list)
save(final_result, file = "/data/zhangh24/multi_ethnic/result/GLGC/CTSLEB/ct_sleb.rdata")

#final_result = read.csv("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/PT.csv")
