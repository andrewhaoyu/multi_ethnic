args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

#i = as.numeric(args[[1]])
l = as.numeric(args[[1]])
j = as.numeric(args[[2]])
v = as.numeric(args[[3]])
#j = as.numeric(args[[3]])
#m = as.numeric(args[[3]])
#i_rep = as.numeric(args[[4]])
#i1 = as.numeric(args[[5]])
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)

eth <- c("EUR","AFR","AMR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
trait = trait_vec[l]
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/1KGLD'),showWarnings = FALSE)
#copy the LD reference data to lscratch
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_eur ",temp.dir,"1KGLD"))
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_afr ",temp.dir,"1KGLD"))
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_amr ",temp.dir,"1KGLD"))
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_eas ",temp.dir,"1KGLD"))
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_sas ",temp.dir,"1KGLD"))
system(paste0("cp /data/zhangh24/software/PRScsx/1KGLD/snpinfo_mult_1kg_hm3 ",temp.dir,"1KGLD"))
data.dir = "/data/zhangh24/multi_ethnic/data/GLGC_cleaned/"
system(paste0("cp ",data.dir,"all_eth.bim ",temp.dir,"all_eth.bim"))


n_vec = rep(0, length(eth))

for(i in 1:5){
  sum = as.data.frame(fread(paste0(data.dir,eth[i],"/",trait,".txt"),header=T))
  n_vec[i] = median(sum$N)
  sum.select = sum %>% 
    rename(SNP=rsID) %>% 
    select(SNP,A1,A2,BETA,P) 
  write.table(sum.select,file = paste0(temp.dir,eth[i],"_sumstats.txt"),row.names = F,col.names = T,quote=F)
}


kg.dir = "/data/zhangh24/KGref_MEGA/GRCh37/"
system(paste0("cp ",kg.dir,eth,"/all_chr.bim ",temp.dir,eth,"all_chr.bim"))
phi = c(1E+00,1E-02,1E-04,1E-6)

path_to_ref = paste0(temp.dir,"1KGLD")
path_to_bim = paste0(temp.dir,"all_eth")
path_to_sum = paste0(temp.dir)

out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/GLGC/prs/PRSCSX/",eth[1],"/",trait,"")




system(paste0("python /data/zhangh24/software/PRScsx/PRScsx.py", 
              " --ref_dir=",path_to_ref,
              " --bim_prefix=",path_to_bim,
              " --sst_file=",path_to_sum,"EUR_sumstats.txt,",path_to_sum,"AFR_sumstats.txt,",path_to_sum,"AMR_sumstats.txt,",path_to_sum,"EAS_sumstats.txt,",path_to_sum,"SAS_sumstats.txt",
              " --n_gwas=",n_vec[1],",",n_vec[2],",",n_vec[3],",",n_vec[4],",",n_vec[5],
              " --pop=EUR,AFR,AMR,EAS,SAS",
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
