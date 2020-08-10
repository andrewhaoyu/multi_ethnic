#select SNPs based on best EUR prs, use coefficients of EUR population
#i represent ethnic group
#l represent causal proportion
#m represent sample size
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
j = as.numeric(args[[4]])
library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01)
#n.snp.mat <- matrix(0,length(pthres),4)
#load EUR r2 result
setwd("/data/zhangh24/multi_ethnic/")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
load(paste0(cur.dir,"LD.clump.result.rdata"))
load("/data/zhangh24/multi_ethnic/result/LD_simulation/r2.mat.eur.rdata")
#keep the EUR results with sample size at 100,000 (m = 4)
r2.mat <- as.data.frame(LD.result.list[[2]]) %>% 
  filter(eth.vec=="EUR"&
           m_vec==4&
           l_vec==l)
#get the best performance eur prs p-value threshold
k = which.max(r2.mat$r2.vec)
LD <- as.data.frame(fread(paste0(cur.dir,eth[1],"/LD_clump_rho_",l,"_size_",4,".clumped")))
clump.snp <- LD[,3,drop=F]  
sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_new/",eth[1],"/summary_out_rho_",l,"_size_",m)))  
colnames(sum.data)[2] <- "SNP"
prs.all <- left_join(clump.snp,sum.data,by="SNP") 
prs.file <- prs.all %>% filter(P<=pthres[k]) %>% 
  select(SNP,A1,BETA)
write.table(prs.file,file = paste0(cur.dir,eth[i],"/prs/prs_eursnp_eurcoef_rho_",l,"_size_",m),col.names = T,row.names = F,quote=F)
system(paste0("/data/zhangh24/software/plink2 --score ",cur.dir,eth[i],"/prs/prs_eursnp_eurcoef_rho_",l,"_size_",m," no-sum no-mean-imputation --bfile ",cur.dir,eth[i],"/chr",j,".tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out ",cur.dir,eth[i],"/prs/prs_eursnp_eurcoef_rho_",l,"_size_",m))

