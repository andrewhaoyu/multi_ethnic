args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#v represent the tuning parameter

j = as.numeric(args[[1]])
v = as.numeric(args[[2]])
#j = as.numeric(args[[3]])
#m = as.numeric(args[[3]])
#i_rep = as.numeric(args[[4]])
#i1 = as.numeric(args[[5]])
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)

eth <- c("EUR","AFR","AMR","EAS","SAS")
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/TEL_dat/"

sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/1KGLD'),showWarnings = FALSE)
#copy the LD reference data to lscratch
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_eur ",temp.dir,"1KGLD"))
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_afr ",temp.dir,"1KGLD"))
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_amr ",temp.dir,"1KGLD"))
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_eas ",temp.dir,"1KGLD"))
system(paste0("cp /data/zhangh24/software/PRScsx/1KGLD/snpinfo_mult_1kg_hm3 ",temp.dir,"1KGLD"))
system(paste0("cp ",data.dir,"all_eth.bim ",temp.dir,"all_eth.bim"))
#out.dir = paste0("/data/zhangh24/multi_ethnic/result/cleaned/clumping_result/",method,"/",eth[i],"/",trait[l],"/")




n.vec = rep(0,4)
#prepare summary statistics
for(i in 1:4){
  load(paste0(data.dir,eth[i],"_sum_data_cleaned.rdata"))
  sum.com.select = sum.com.select %>% 
    filter(CHR==j)
  sum.select = sum.com.select %>% 
    select(SNP, A1, A2, BETA, P) 
  n.vec[i] = median(sum.com.select$N)
  write.table(sum.select,file = paste0(temp.dir,eth[i],"_sumstats.txt"),row.names = F,col.names = T,quote=F)
}

phi = c(1E+00,1E-02,1E-04,1E-6)

path_to_ref = paste0(temp.dir,"1KGLD")
path_to_bim = paste0(temp.dir,"all_eth")
path_to_sum = paste0(temp.dir)

out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/TEL/prscsx")




system(paste0("python /data/zhangh24/software/PRScsx/PRScsx.py", 
              " --ref_dir=",path_to_ref,
              " --bim_prefix=",path_to_bim,
              " --sst_file=",path_to_sum,"EUR_sumstats.txt,",path_to_sum,"AFR_sumstats.txt,",path_to_sum,"AMR_sumstats.txt,",path_to_sum,"EAS_sumstats.txt",
              " --n_gwas=",n.vec[1], ",", n.vec[2], ",", n.vec[3], ",", n.vec[4],
              " --pop=EUR,AFR,AMR,EAS",
              " --chrom=",j,
              " --phi=",phi[v],  
              " --out_dir=",out.dir.prs,
              " --out_name=sum"))




# out.file = paste0("/data/zhangh24/multi_ethnic/result/cleaned/organize_prs/PRSCSx/")
# for(i in 1:5){
#   system(paste0("mkdir ",out.file,eth[i],"/"))
#   for(l in 1:7){
#     system(paste0("mkdir ",out.file,eth[i],"/",trait[l]))  
#   }
# }
