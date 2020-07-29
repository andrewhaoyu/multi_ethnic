#goal extract causal SNPs
#since GCTA has limited read-in speed compared to plink
#first extract the causal SNPs
#then simulate phenotypes
args <- commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/cau.snp.infor.list.rdata")
eth <- c("EUR","AFR","AMR","EAS","SAS")
cau.snp.infor <- cau.snp.infor.list[[1]]
library(data.table)
EUR.bi <- cau.snp.infor$EUR>=0.01&
  cau.snp.infor$EUR<=0.99
AFR.bi <-  cau.snp.infor$AFR>=0.01&
  cau.snp.infor$AFR<=0.99
AMR.bi <-  cau.snp.infor$AMR>=0.01&
  cau.snp.infor$AMR<=0.99
EAS.bi <-  cau.snp.infor$EAS>=0.01&
  cau.snp.infor$EAS<=0.99
SAS.bi <-  cau.snp.infor$SAS>=0.01&
  cau.snp.infor$SAS<=0.99
eth.bi <- cbind(EUR.bi,AFR.bi,AMR.bi,EAS.bi,SAS.bi)
library(dplyr)
MAF <- cau.snp.infor %>% select(EUR,AFR,AMR,EAS,SAS)
MAF <- as.data.frame(MAF)
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
  idx <- which(eth.bi[,i]==1)
  select.cau.snp <- cbind(cau.snp.infor[idx,1])
  write.table(select.cau.snp,file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i],"/select.cau.snp.txt"),row.names = F,col.names = F,quote=F)

  system(paste0("/data/zhangh24/software/plink2 --bfile ",cur.dir,eth[i],"/all_chr.tag --extract ",cur.dir,eth[i],"/select.cau.snp.txt --make-bed --out select.cau.snp"))
