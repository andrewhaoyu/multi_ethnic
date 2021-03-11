#LD_clumping for ARIC data
#load the p-value results
args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

#i = as.numeric(args[[1]])
#l = as.numeric(args[[2]])

j = 1
i = 1

l = 1
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)

eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")
library(bigsnpr)
# for(i in 1:2){
#   for(l in c(1,3)){
setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic")
temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
load(paste0(temp.dir,"all_keep_chr_",j,".rdata"))
load(paste0(temp.dir,"multi_PRS_chr_",j,".rdata"))
