load( "/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/prediction_summary.rdata")

temp = prediction.result %>% filter(method_vec%in%c("CT-SLEB (two ancestries)"))
temp2 = prediction.result %>% filter(method_vec%in%c("CT-SLEB (five ancestries)"))
temp-temp2


temp$new_reuslt = temp2$result-temp$result
idx <- which(temp$new_reuslt<0)
temp[idx,]
