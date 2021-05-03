dit = 3
result.table.list = list()
library(data.table)
l = 1

for(l in 1:3){
  load(paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/result/auc_result_",l))
  result.table = data.frame()
  for(k in 1:3){
    result.table[1,k] = paste0(round(auc.result[[1]][k],dit),
                             " (",
                             round(auc.result[[2]][k],dit),
                             ")")
                             
  }
  result.table.list[[l]] = result.table
}



result.table = rbindlist(result.table.list)
colnames(result.table) = names(auc.result[[1]])
write.csv(result.table,file = "/data/zhangh24/multi_ethnic/result/breast_cancer/result/auc_result_best_eur.csv")
