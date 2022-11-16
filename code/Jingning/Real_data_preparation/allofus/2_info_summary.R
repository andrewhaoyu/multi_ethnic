
#trait='HDL'
#ethnic='AMR'


library(readr)
library(bigreadr)
library(dplyr)


eth_all <- c("AFR","AMR","EUR")
trait_all <- c("height","bmi")

res <- tibble()
for (t in 1:length(trait_all)){
  trait <- trait_all[t]

res_t <- matrix(nrow=length(eth_all),ncol=4)

for (i in 1:length(eth_all)){
  ethnic <- eth_all[i]

dat <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/allofus_cleaned/",ethnic,"/",trait,".txt"))
res_t[i,1] <- max(dat$N)

num_rec <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/allofus_cleaned/",ethnic,"/",trait,"_NofSNP.rds"))

  res_t[i,2] <- num_rec["original"]
  res_t[i,3] <- num_rec["mega"]
  #res_t[i,4] <- num_rec["remove_duplicated_rsid"] # remove missing p-value and uplicated rsid
  #res_t[i,5] <- num_rec["remove_small_sample_size"] # remove remove_anbiguous_strand and small sample size
  #res_t[i,6] <- num_rec["in1000G"]
  res_t[i,4] <- num_rec["common1%"]

  }

  res_t <- data.frame(ethnic=eth_all,trait=trait, res_t)
  colnames(res_t) <- c("ethnic","trait","sample size","Total Number of SNPs","SNPs overlapping with HM3 + MEGA chips","SNPs overlapping with HM3 + MEGA chips and  passing quality control")

  res <- rbind(res, res_t)
}

write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/summary_data/allofus_cleaned/info_summary.txt"))



