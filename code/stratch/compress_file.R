#goal: compress file
eth_vec <- c("EUR","AFR","AMR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
setwd("/data/zhangh24/multi_ethnic/result/GLGC/pgs")
library(data.table)
for(i in c(1,2,4,5)){
  for(l in 1:4){
    eth = eth_vec[i]
    trait = trait_vec[l]
    data = fread(file = paste0(trait,"_",eth,"_","LDpred2.txt"))
    out_filename = paste0("/data/zhangh24/multi_ethnic/result/GLGC/pgs_catalog/",
                          trait,"_",eth,"_","LDpred2.txt.gz")
    write.table(data, file = gzfile(out_filename), sep = "\t", 
                row.names = FALSE, quote = FALSE, col.names = TRUE)

  }
}


for(i in c(2,4,5)){
  for(l in 1:4){
    eth = eth_vec[i]
    trait = trait_vec[l]
    data = fread(file = paste0(trait,"_",eth,"_","weighted_LDpred2.txt"))
    out_filename = paste0("/data/zhangh24/multi_ethnic/result/GLGC/pgs_catalog/",
                          trait,"_",eth,"_","weighted_LDpred2.txt.gz")
    write.table(data, file = gzfile(out_filename), sep = "\t", 
                row.names = FALSE, quote = FALSE, col.names = TRUE)
    
  }
}



eth_vec <- c("EUR","AFR","AMR")
trait_vec <-c("height","bmi")

setwd("/data/zhangh24/multi_ethnic/result/AOU/pgs")

for(i in c(1:2)){
  for(l in 1:2){
    eth = eth_vec[i]
    trait = trait_vec[l]
    data = fread(file = paste0(trait,"_",eth,"_","LDpred2.txt"))
    out_filename = paste0("/data/zhangh24/multi_ethnic/result/AOU/pgs_catalog/",
                          trait,"_",eth,"_","LDpred2.txt.gz")
    write.table(data, file = gzfile(out_filename), sep = "\t", 
                row.names = FALSE, quote = FALSE, col.names = TRUE)
    
  }
}


for(i in c(2)){
  for(l in 1:2){
    eth = eth_vec[i]
    trait = trait_vec[l]
    data = fread(file = paste0(trait,"_",eth,"_","weighted_LDpred2.txt"))
    out_filename = paste0("/data/zhangh24/multi_ethnic/result/AOU/pgs_catalog/",
                          trait,"_",eth,"_","weighted_LDpred2.txt.gz")
    write.table(data, file = gzfile(out_filename), sep = "\t", 
                row.names = FALSE, quote = FALSE, col.names = TRUE)
    
  }
}
