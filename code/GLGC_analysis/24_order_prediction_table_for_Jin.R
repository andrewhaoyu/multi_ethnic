load("/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/GLGC/analysis_result/prediction.result.summary.rdata")
order.trait =  c("HDL",
                 "LDL",
                 "logTG",
                 "TC")
library(dplyr)
load("/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/GLGC/analysis_result/weighted_prs_five_ans.rdata")
head(final_result)
final_result_update = final_result %>% 
  rename(result = r2) %>% 
  mutate(method_vec = method)
prediction.result.com = rbind(prediction.result,final_result_update)

load("/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/AOU/analysis_result/prediction.result.summary.rdata")

load("/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/AOU/analysis_result/weighted_prs_three_ans.rdata")
head(final_result)
final_result_update = final_result %>% 
  rename(result = r2) %>% 
  mutate(method_vec = method)
prediction.result.com2 = rbind(prediction.result,final_result_update)

prediction.result.all = rbind(prediction.result.com,
                              
                              prediction.result.com2)
idx <- which(duplicated(prediction.result.all))
save(prediction.result.all,
    file =  "/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/GLGC/analysis_result/prediction.result.summary.combined.rdata")
write.csv(prediction.result.all, file = "/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/GLGC/analysis_result/prediction.result.summary.combined.csv",
          quote = F)
idx = which(duplicated(prediction.result.all))
