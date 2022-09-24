#goal: preprocess aou summary statistics
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])

library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)

eth_vec <- c("EUR","AFR","AMR")
trait_vec <-c("height","bmi")
eth = eth_vec[i]
trait = trait_vec[l]
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/AOU_cleaned/"
out.dir = paste0("/data/zhangh24/multi_ethnic/result/AOU/clumping_result/PT/",eth,"/",trait,"/")
#load grch37 and 38 match information
snp_infor = readRDS("/data/zhangh24/multi_ethnic/result/LD_simulation_new/SNP_GRCh37_38_match_update.rds")
snp_infor_select = snp_infor %>% 
  select(rsid, pos37) %>% 
  rename(rsID = rsid)
#########LD clumping#################
#load gwas summary statistics
sum.data = as.data.frame(fread(paste0(data.dir,eth,"/",trait,".txt"),header=T))
sum_data_update = left_join(sum.data, snp_infor_select)
sum_data_update = sum_data_update %>% na.omit()
write.table(sum_data_update, file = paste0(data.dir,eth,"/",trait,"_update.txt"), 
            row.names = F, col.names = T, quote = F)
