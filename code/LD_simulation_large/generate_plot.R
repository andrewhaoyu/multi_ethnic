setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_new")
library(ggplot2)
load("LD.clump.result.rdata")
sample_size <- rep(c("15000","45000","80000"),15)
LD.clump.result <- cbind(LD.clump.result,sample_size)

p_list <- list()