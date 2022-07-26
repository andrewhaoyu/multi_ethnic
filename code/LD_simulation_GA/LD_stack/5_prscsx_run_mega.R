#use PRS-csx on simulation
args = commandArgs(trailingOnly = T)
i = 2
l = 1
m = 1
i_rep = as.numeric(args[[1]])
i1 = 1
j = as.numeric(args[[2]])
k = as.numeric(args[[3]])
library(data.table)
library(tidyverse)
eth <- c("EUR","AFR","AMR","EAS","SAS")
sid<-Sys.getenv('SLURM_JOB_ID')
#system(paste0("rm -rf /lscratch/",sid,'/test'))
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/1KGLD'),showWarnings = FALSE)
dir.create(paste0('/lscratch/',sid,'/test/1KGLD/ldblk_1kg_eur'),showWarnings = FALSE)
dir.create(paste0('/lscratch/',sid,'/test/1KGLD/ldblk_1kg_',tolower(eth[i])),showWarnings = FALSE)

#copy the LD reference data to lscratch
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD_MEGA/ldblk_1kg_eur/ldblk_1kg_chr",j,".hdf5 ",temp.dir,"1KGLD/ldblk_1kg_eur/"))
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD_MEGA/ldblk_1kg_",tolower(eth[i]),"/ldblk_1kg_chr",j,".hdf5 ",temp.dir,"1KGLD/ldblk_1kg_",tolower(eth[i]),"/"))
system(paste0("cp /data/zhangh24/software/PRScsx/1KGLD_MEGA/snpinfo_mult_1kg_hm3 ",temp.dir,"1KGLD"))
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
prs_cs_ref = as.data.frame(fread("/data/zhangh24/software/PRScsx/1KGLD_MEGA/snpinfo_mult_1kg_hm3",header=T))
prs_cs_ref = prs_cs_ref %>% select(SNP)
#match the summary stat with prs csx reference
summary.eur = inner_join(summary.eur,prs_cs_ref,by = c("rs_id"="SNP"))
#A1 is the effect allele
#A2 is the non-effect allele; generate A2 for prs-csx
summary.eur.sub = summary.eur %>% 
  separate(SNP,into = c("non_im","non_im2","first_allele","second_allele"),
           remove = F) %>% 
  mutate(A2 =ifelse(A1==second_allele,first_allele,second_allele)) 
#filter the SNP to CHR j
summary.eur.select = summary.eur.sub %>% 
  filter(CHR == j) %>% 
  #select(rs_id,A1,A2,BETA,P,CHR) %>% 
  select(rs_id, A1, A2, BETA, P) %>% 
  rename(SNP = rs_id)
write.table(summary.eur.select,file = paste0(temp.dir,"EUR_sumstats.txt"),row.names = F,col.names = T,quote=F)
#merge with prs-csx reference data
#prepare target summary stat for a specific chr
summary.tar <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  

summary.tar = left_join(summary.tar,snp.infor,by="SNP")
#merge with prs-csx reference data
summary.tar = inner_join(summary.tar,prs_cs_ref,by = c("rs_id"="SNP"))
#A1 is the effect allele
#A2 is the non-effect allele; generate A2 for prs-csx
summary.tar.sub = summary.tar %>% 
  separate(SNP,into=c("non_im","non_im2","first_allele","second_allele"),
           remove=F) %>% 
  mutate(A2 =ifelse(A1==second_allele,first_allele,second_allele)) 

#filter the SNP to CHR j
summary.tar.select = summary.tar.sub %>% 
  filter(CHR == j) %>% 
  select(rs_id,A1,A2,BETA,P) %>% 
  #select(rs_id,A1,A2,BETA,P,CHR) %>% 
  rename(SNP=rs_id)
write.table(summary.tar.select,file = paste0(temp.dir,eth[i],"_sumstats.txt"),row.names = F,col.names = T,quote=F)
bim = read.table(paste0(cur.dir,eth[i],"/all_chr_test.mega.bim"))
bim.update =inner_join(summary.tar,bim,by = c("SNP"="V2")) %>% 
  select(V1,rs_id,V3,V4,V5,V6)
bim.select = bim.update %>% 
  filter(V1 == j)

write.table(bim.update,file = paste0(temp.dir,"all_chr_test.mega.bim"),row.names = F,col.names = F,quote=F)
#phi = c(1E-6,1E-4)
#phi = c(1E-02)
phi = c(1,1E-02,1E-04,1E-06)
# for(k in 1:length(phi)){
#   for(j in 1:22){
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
                  " --out_dir=",out.dir.sum,eth[i],"/prscsx_mega",
                  " --out_name=rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
    

    # summary.eur.select = summary.eur.sub %>% 
    #   filter(CHR==1) %>% 
    #   select(SNP,A1,A2,BETA,P)
    # 
    # summary.eur.select = summary.eur.select[50272:50272,]
    # 
    # write.table(summary.eur.select,file = paste0(temp.dir,"EUR_sumstats.txt"),row.names = F,col.names = T,quote=F)
    # 
    # 
    # summary.tar.select = summary.tar.sub %>% 
    #   filter(CHR==1) %>% 
    #   select(SNP,A1,A2,BETA,P)
    # summary.tar.select = summary.tar.select[50272:50272,]
    # 
    # write.table(summary.tar.select,file = paste0(temp.dir,eth[i],"_sumstats.txt"),row.names = F,col.names = T,quote=F)
    # 
    # 
    # system(paste0("export MKL_NUM_THREADS=2; export NUMEXPR_NUM_THREADS=2; export OMP_NUM_THREADS=2;
    #             python /data/zhangh24/software/PRScsx/PRScsx.py", 
    #               " --ref_dir=",path_to_ref,
    #               " --bim_prefix=",path_to_bim,
    #               " --sst_file=",path_to_sum,"EUR_sumstats.txt,",path_to_sum,eth[i],"_sumstats.txt",
    #               " --n_gwas=100000,",size_list[m],
    #               " --pop=EUR,",eth[i],
    #               " --chrom=",j,
    #               " --phi=",phi[k], 
    #               " --out_dir=",out.dir.sum,eth[i],"/prscsx_mega",
    #               " --out_name=rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
    # 
    # write.table(summary.eur.select,file = paste0("/data/zhangh24/test1/PRScsx/test_data/","EUR_sumstats.txt"),row.names = F,col.names = T,quote=F)
    # write.table(summary.tar.select,file = paste0("/data/zhangh24/test1/PRScsx/test_data/",eth[i],"_sumstats.txt"),row.names = F,col.names = T,quote=F)
    # 
    # bim.update.select = bim.update %>% 
    #   filter(V1==1)
    # write.table(bim.update.select,file = paste0("/data/zhangh24/test1/PRScsx/test_data/","test.bim"),row.names = F,col.names = F,quote=F)
    # 
    # 
    # system(paste0("cd /data/zhangh24/test1/; ",
    #               "python ./PRScsx/PRScsx.py --ref_dir=./ref2 --bim_prefix=./PRScsx/test_data/test --sst_file=./PRScsx/test_data/EUR_sumstats.txt,./PRScsx/test_data/AFR_sumstats.txt --n_gwas=100000,15000 --pop=EUR,AFR --chrom=1 --phi=1e-02 --out_dir=./ --out_name=test2"))
    # 
    # 
    # 