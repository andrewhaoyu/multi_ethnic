eth_vec <- c("EUR","AFR","AMR")
trait_vec <-c("height","bmi")

result = list()
temp = 1
for(l in 1:2){
  for(i in 1:3){
    eth = eth_vec[i]
    trait = trait_vec[l]
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/AOU/clumping_result/PT/",eth,"/",trait,"/")
    load(paste0(out.dir, "CT.result"))
    result[[temp]] = ct.result[[1]]
    temp = temp+ 1
  }
}

final_result = rbindlist(result)
save(final_result, file = "/data/zhangh24/multi_ethnic/result/AOU/clumping_result/PT/PT.rdata")

#final_result = read.csv("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/PT.csv")
