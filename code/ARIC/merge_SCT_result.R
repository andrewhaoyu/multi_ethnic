eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")
library(bigsnpr)
# for(i in 1:2){
#   for(l in c(1,3)){
setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic")
total = 4
trait_vec = rep("c",total)
eth_vec = rep("c",total)
r2 = rep(0,total)
temp = 1
for(l in c(1,3)){
  for(i in 1:2){
    temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
    data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
    geno.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[1],"/",eth[i],"/")
    out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
    load(paste0(out.dir,"SCT_r2.rdata"))
    eth_vec[temp] = eth[i]
    trait_vec[temp] = trait[l]
    r2[temp] = r2.list[[1]]
    temp = temp+1
  }
  
}
ARIC.result.SCT = data.frame(eth_vec,trait_vec,
                             r2,method_vec = "SCT")
save(ARIC.result.SCT,file = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/ARIC.result.SCT.rdata"))
  