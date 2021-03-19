eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")
library(bigsnpr)
# for(i in 1:2){
#   for(l in c(1,3)){
setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic")



for(l in c(1,3)){
  for(i in 1:2){
    temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
    data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
    out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
    j = 1
    load(paste0(temp.dir,"multi_PRS_chr_",j,".rdata.gzip"))
    multi_PRS_final = multi_PRS
    for(j in 2:22){
      load(paste0(temp.dir,"multi_PRS_chr_",j,".rdata.gzip")) 
      multi_PRS_final = multi_PRS_final + multi_PRS
    }
    multi_PRS = multi_PRS_final
    save(multi_PRS, file = paste0(temp.dir,"multi_PRS.rdata"))
  }
}





#paste0(temp.dir,"multi_PRS_chr_",j,".rdata.gzip")