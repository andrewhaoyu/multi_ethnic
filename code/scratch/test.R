
library(data.table)
snp_infor = readRDS("/gpfs/gsfs12/users/BB_Bioinformatics/ProjectData/SNP_infor_GRCh37_38/SNP_GRCh37_38_match_update.rds")
library(dplyr)
#load smkces data
smkces = read.table(gzfile("/gpfs/gsfs11/users/zhangh24/GSCAN/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_MULTI.txt.gz"),header= T,sep = "\t")
#load mds data
mds = read.table(gzfile("/gpfs/gsfs11/users/zhangh24/GSCAN/GSCAN_SmkCes_2022_MDS_MULTI.txt.gz"),header=T)
#find SNP rs2459985 
idx = which(smkces$POS==2067849)
smkces[idx,]
#add intercept
mds = cbind(1,mds)
gamma_vec = c(-0.013219, -0.000166,5.4e-05,
              0.000227,0.000953)
se_vec = c(0.003341,5.8e-05,0.000154,0.000222,0.000867)
#calcualte beta
beta_vec = as.matrix(mds)%*%gamma_vec
#fixed effect-meta analyses
mds2 = mds^2
var_vec = as.matrix(mds2)%*%se_vec^2

w_vec = 1/var_vec
sum(beta_vec*w_vec)/(sum(w_vec))
