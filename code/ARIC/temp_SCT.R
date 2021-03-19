eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")
library(bigsnpr)
# for(i in 1:2){
#   for(l in c(1,3)){
setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/dcl02/leased/chatterj/hzhang1/multi_ethnic/result/ARIC/eGFRcr/EUR/"
#out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
snp_readBed(paste0(data.dir,"all_chr.bed"))
