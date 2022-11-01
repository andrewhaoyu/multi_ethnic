#use PRS-csx on simulation
args = commandArgs(trailingOnly = T)
l = as.numeric(args[[1]])
m = as.numeric(args[[2]])
i_rep = as.numeric(args[[3]])
i1 = as.numeric(args[[4]])
k = as.numeric(args[[5]])
j = as.numeric(args[[6]])
library(data.table)
library(tidyverse)
eth <- c("EUR","AFR","AMR","EAS","SAS")
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
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
system(paste0("cp ",out.dir, "all_eth.bim ", temp.dir,"all_eth.bim"))
path_to_bim = paste0(temp.dir,"all_eth")
setwd("/data/zhangh24/multi_ethnic/")
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
snp.infor = snp.infor.match %>% 
  rename(SNP=id) %>% 
  select(SNP,rs_id)

prs_cs_ref = read.table("/gpfs/gsfs11/users/zhangh24/software/PRScsx/1KGLD/snpinfo_mult_1kg_hm3",header=T)
prs_cs_ref = prs_cs_ref %>% select(SNP)

#prepare summary statistics
for(i in 1:5){
  if(i==1){
    sum = as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))    
  }else{
    sum = as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
  }
  sum.update = left_join(sum,snp.infor,by="SNP")
  sum.update = inner_join(sum.update,prs_cs_ref,by = c("rs_id"="SNP"))
  sum.update = sum.update %>% 
    separate(SNP,into=c("non_im","non_im2","first_allele","second_allele"),
             remove=F) %>% 
    mutate(A2 =ifelse(A1==second_allele,first_allele,second_allele)) %>% 
    select(rs_id,A1,A2,BETA,P) %>% 
    rename(SNP=rs_id)
  
  write.table(sum.update,file = paste0(temp.dir,eth[i],"_sumstats.txt"),row.names = F,col.names = T,quote=F)
}


#phi = c(1E-6,1E-4)
#phi = c(1E-02)
phi = c(1,1E-02,1E-04,1E-06)
#for(k in 1:length(phi)){
#for(j in 1:22){
  print(j)
  #create bim file for prscsx
  #run prs-csx
  path_to_ref = paste0(temp.dir,"1KGLD")
  path_to_bim = paste0(temp.dir,"all_eth")
  path_to_sum = paste0(temp.dir)
  size_list = c("15000","45000","80000","100000")
  
  system(paste0("export MKL_NUM_THREADS=2; export NUMEXPR_NUM_THREADS=2; export OMP_NUM_THREADS=2;
                python /data/zhangh24/software/PRScsx/PRScsx.py", 
                " --ref_dir=",path_to_ref,
                " --bim_prefix=",path_to_bim,
                " --sst_file=",path_to_sum,"EUR_sumstats.txt,",path_to_sum,"AFR_sumstats.txt,",path_to_sum,"AMR_sumstats.txt,",path_to_sum,"EAS_sumstats.txt,",path_to_sum,"SAS_sumstats.txt",
                " --n_gwas=100000,",size_list[m],",",size_list[m],",",size_list[m],",",size_list[m],
                " --pop=EUR,AFR,AMR,EAS,SAS",
                " --chrom=",j,
                " --phi=",phi[k], 
                " --out_dir=",out.dir.sum,eth[1],"/prscsx",
                " --out_name=update_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
  
  
  
#}

#}


