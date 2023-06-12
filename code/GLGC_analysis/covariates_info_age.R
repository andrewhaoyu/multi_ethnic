setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/GLGC_cleaned/"
out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/",eth,"/",trait,"/")
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)
pheno.dir = "/data/zhangh24/multi_ethnic/data/UKBB/phenotype/"
eth_vec <- c("AFR","AMR","EAS","EUR","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
covar_list = list()
l = 1
trait = trait_vec[l]
result = rep("c",5)
result_sex = rep(0,5)
for(i in 1:5){
  eth = eth_vec[i]
  covar <- as.data.frame(fread(paste0(pheno.dir,"/covariates/tuning+validation/",eth,"_all_data.txt")))  
  colnames(covar) = c('id','sex','age',paste0('pc',1:10))
  result[i] = paste0(round(mean(covar$age),2)," (",
                           round(sd(covar$age),2),")")
  result_sex[i] = round(sum(1-covar$sex)/nrow(covar),2)
  covar_list[[i]] = covar
}

all_covar = rbindlist(covar_list)
result = c(paste0(round(mean(all_covar$age),2)," (",
             round(sd(all_covar$age),2),")"),result)
result_sex = c(paste0(round(mean(covar$sex),2)),result_sex)
print(result)
result_sex
