eth_vec <- c("EUR","AFR","AMR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
library(data.table)
result = list()
temp = 1
  for(l in 1:4){
    for(i in 1:5){
      eth = eth_vec[i]
      trait = trait_vec[l]
      out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/",eth,"/",trait,"/")
      load(paste0(out.dir, "CT.result"))
      result[[temp]] = ct.result[[1]]
      temp = temp+ 1
  }
  }

final_result = rbindlist(result)
save(final_result, file = "/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/PT.rdata")

#final_result = read.csv("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/PT.csv")
