#calculate the AUC for EUR SNPs with EUR coefficients
#with target ethnic coefficients
#with EB coefficients
#calculate AUC for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
#q for three different methods

#i_rep = as.numeric(args[[5]])


library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
#n <- 120000
n.train.vec <- c(15000,45000,80000,100000)
n.rep = 10

#r2 mat represent the r2 matrix for the testing dataset
#column represent the ethnic groups
#row represent different p-value threshold

cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
setwd("/data/zhangh24/multi_ethnic/")


result.list = list()
temp = 1
for(i in 2:5){
  for(i1 in 1:5){
    for(l in 1:3){
      load(paste0(out.dir,
                  "prscsx_five.result_sub_ga_",i1,"_eth_",i,"_rho_",l,".rdata"))
      result.list[[temp]] = prscsx.result
      temp = temp + 1
#load result
# }
    }


  }
}
prscsx.result = rbindlist(result.list)
save(prscsx.result,file = paste0(out.dir,
                                 "prscsx_five.result.rdata"))




