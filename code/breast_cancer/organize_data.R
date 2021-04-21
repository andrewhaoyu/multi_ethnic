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

library(dplyr)
library(data.table)
load("/dcs04/nilanjan/data/23andme/snpinfo/all_snp_info.RData")
snpinfo.temp = snpinfo
# idx <- which(is.na(snpinfo$im.data.id))
# snpinfo[idx[1:10],]
mega <- fread("/dcs04/nilanjan/data/23andme/cleaned/mega+hapmap3_autosomal_in_1000G",header=F)
mega.match <- inner_join(mega,snpinfo,by = c("V1"="assay.name"))
idx <- which(is.na(mega.match$im.data.id))
print(length(idx))


load("/dcs04/nilanjan/data/23andme/cleaned/snpinfo/snpinfo_mega.RData")
idx <- which(mega$V1%in%snpinfo_mega$assay.name==F)

