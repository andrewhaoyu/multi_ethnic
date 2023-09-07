eth_vec <- c("EUR","AFR","AMR")
trait_vec <-c("height","bmi")

result_list = list()
temp = 1
for(l in 1:2){
  for(i in 1:3){
    eth = eth_vec[i]
    trait = trait_vec[l]
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/AOU/weighted_prs/",eth,"/",trait,"/")
    load(paste0(out.dir, "weighted_prs_ct_three_ans.result"))
    result_list[[temp]] = result %>% mutate(method = "Weighted PRS (C + T) three ancestries")
    temp = temp+ 1
  }
}

final_result = rbindlist(result_list)
save(final_result, file = "/data/zhangh24/multi_ethnic/result/AOU/weighted_prs/weighted_prs_three_ans.rdata")

#final_result = read.csv("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/PT.csv")
