#organize the summarize data for breast cancer
setwd("/data/zhangh24/multi_ethnic/data/AABC_data")
library(data.table)
trait = c("overall","ERpos","ERneg")
set = c("train","valid")

for(k in 1:3){
  for(l in 1:2){
    
    data.list = list()
    for(j in 1:22){
      print(c(k,l,j))
      filename = paste0("RAGOC.MAF005R7.",trait[k],".",set[l],".chr",j,".logistic")
      sum.data.sub <- as.data.frame(fread(filename))
      data.list[[j]] = sum.data.sub
    }  
    sum.data = rbindlist(data.list)
    save(sum.data,file = paste0("BC_AFR_",trait[k],"_",set[l],".rdata"))
  }
  
}

