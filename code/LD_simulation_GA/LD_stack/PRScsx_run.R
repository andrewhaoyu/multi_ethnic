#use PRS-csx on simulation
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])
i1 = as.numeric(args[[5]])
library(data.table)
library(tidyverse)
eth <- c("EUR","AFR","AMR","EAS","SAS")
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/1KGLD'),showWarnings = FALSE)
#copy the LD reference data to lscratch
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_eur ",temp.dir,"1KGLD"))
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_",tolower(eth[i])," ",temp.dir,"1KGLD"))
system(paste0("cp /data/zhangh24/software/PRScsx/1KGLD/snpinfo_mult_1kg_hm3 ",temp.dir,"1KGLD"))
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
summary.eur <- as.data.frame(fread(paste0(out.dir.sum,eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))  
#prepare EUR summary stat for a specific chr
#load 1kg snp list
setwd("/data/zhangh24/multi_ethnic/")
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
snp.infor = snp.infor.match %>% 
  rename(SNP=id) %>% 
  select(SNP,rs_id)
summary.eur = left_join(summary.eur,snp.infor,by="SNP")
prs_cs_ref = read.table("/gpfs/gsfs11/users/zhangh24/software/PRScsx/1KGLD/snpinfo_mult_1kg_hm3",header=T)
prs_cs_ref = prs_cs_ref %>% select(SNP)
summary.eur = inner_join(summary.eur,prs_cs_ref,by = c("rs_id"="SNP"))
summary.eur = summary.eur %>% 
  separate(SNP,into=c("non_im","non_im2","first_allele","second_allele"),
           remove=F) %>% 
  mutate(A2 =ifelse(A1==second_allele,first_allele,second_allele)) %>% 
  select(rs_id,A1,A2,BETA,P) %>% 
  rename(SNP=rs_id)
write.table(summary.eur,file = paste0(temp.dir,"EUR_sumstats.txt"),row.names = F,col.names = T,quote=F)
#merge with prs-csx reference data
#prepare target summary stat for a specific chr
summary.tar <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
summary.tar = left_join(summary.tar,snp.infor,by="SNP")
#merge with prs-csx reference data
summary.tar = inner_join(summary.tar,prs_cs_ref,by = c("rs_id"="SNP"))
summary.tar.sub = summary.tar %>% 
  separate(SNP,into=c("non_im","non_im2","first_allele","second_allele"),
           remove=F) %>% 
  mutate(A2 =ifelse(A1==second_allele,first_allele,second_allele)) %>% 
  select(rs_id,A1,A2,BETA,P) %>% 
  rename(SNP=rs_id)
write.table(summary.tar.sub,file = paste0(temp.dir,eth[i],"_sumstats.txt"),row.names = F,col.names = T,quote=F)
bim = read.table(paste0(cur.dir,eth[i],"/all_chr_test.mega.bim"))
bim.update =inner_join(summary.tar,bim,by = c("SNP"="V2")) %>% 
  select(V1,rs_id,V3,V4,V5,V6)
write.table(bim.update,file = paste0(temp.dir,"all_chr_test.mega.bim"),row.names = F,col.names = F,quote=F)
#phi = c(1E-6,1E-4)
#phi = c(1E-02)
phi = c(1,1E-02,1E-04,1E-06)
for(k in 1:length(phi)){
  for(j in 1:22){
    print(j)
    #create bim file for prscsx
    #run prs-csx
    path_to_ref = paste0(temp.dir,"1KGLD")
    path_to_bim = paste0(temp.dir,"all_chr_test.mega")
    path_to_sum = paste0(temp.dir)
    size_list = c("15000","45000","80000","100000")
    system(paste0("export MKL_NUM_THREADS=2; export NUMEXPR_NUM_THREADS=2; export OMP_NUM_THREADS=2;
                python /data/zhangh24/software/PRScsx/PRScsx.py", 
                  " --ref_dir=",path_to_ref,
                  " --bim_prefix=",path_to_bim,
                  " --sst_file=",path_to_sum,"EUR_sumstats.txt,",path_to_sum,eth[i],"_sumstats.txt",
                  " --n_gwas=100000,",size_list[m],
                  " --pop=EUR,",eth[i],
                  " --chrom=",j,
                  " --phi=",phi[k], 
                  " --out_dir=",out.dir.sum,eth[i],"/prscsx",
                  " --out_name=rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
    
    
    
  }
  
}


