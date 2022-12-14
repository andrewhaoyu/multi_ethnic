eth_vec <- c("EUR","AFR","AMR","EAS", "SAS")
trait_vec <-c("HDL","LDL",
              "logTG",
              "TC")

result_list = list()
temp = 1
for(l in 1:4){
  for(i in 1:5){
    eth = eth_vec[i]
    trait = trait_vec[l]
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/weighted_prs/",eth,"/",trait,"/")
    load(paste0(out.dir, "weighted_prs_ct_five_ans.result"))
    result_list[[temp]] = result
    temp = temp+ 1
  }
}

final_result = rbindlist(result_list)
save(final_result, file = "/data/zhangh24/multi_ethnic/result/GLGC/weighted_prs/weighted_prs_five_ans.rdata")

#final_result = read.csv("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/PT.csv")
