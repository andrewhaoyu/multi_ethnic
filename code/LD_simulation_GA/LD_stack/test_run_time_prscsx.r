args = commandArgs(trailingOnly = T)
t_rep = args[[1]]
k = args[[2]]
#test the running time and memory for TDLD-SLEB and PRS-CSx on CHR 22
time1 = proc.time()
i = 2
phi = c(1,1E-02,1E-04,1E-06)
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
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/time_mem/"
phi = c(1,1E-02,1E-04,1E-06)
m = 1;j = 22
path_to_ref = paste0(temp.dir,"1KGLD")
path_to_bim = paste0(out.dir,"AFR_ref_chr22_rs_id")
path_to_sum = paste0(out.dir)
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
              " --out_dir=",out.dir,
              " --out_name=test"))
#calculate prs for the two coefficients
phi_ch = c("1e+00","1e-02","1e-04","1e-06")
#calcualte prs for EUR
prs = fread(paste0(out.dir,"test_",eth[i],
                   "_pst_eff_a1_b0.5_phi",phi_ch[k],"_chr22.txt"))
prs = prs %>% select(V1,V2,V6)
#load AFR summary data
load(paste0(out.dir,"tar_sumdata.rdata"))
prs.match = left_join(prs,summary.tar,by=c("V2"="rs_id"))
dir.create(paste0('/lscratch/',sid,'/test/prs/'),showWarnings = FALSE)
temp.dir.prs = paste0('/lscratch/',sid,'/test/prs/')
system(paste0("cp ",out.dir,"AFR_ref_chr22.bed ",temp.dir,"AFR_ref_chr22.bed"))
system(paste0("cp ",out.dir,"AFR_ref_chr22.bim ",temp.dir,"AFR_ref_chr22.bim"))
system(paste0("cp ",out.dir,"AFR_ref_chr22.fam ",temp.dir,"AFR_ref_chr22.fam"))
prs.file <- prs.match %>% 
  select(SNP,A1,BETA)
colSums(is.na(prs.file))
write.table(prs.file,file = paste0(temp.dir.prs,"prs_file"),col.names = T,row.names = F,quote=F)

res = system(paste0("/data/zhangh24/software/plink2_alpha --score-col-nums 3 --threads 2 --score ",temp.dir.prs,"prs_file  header no-mean-imputation --bfile ",temp.dir,"AFR_ref_chr22  --out ",temp.dir.prs,"prs_csx_",eth[i],"_phi",phi[k]))
prs = fread(paste0(out.dir,"test_",eth[1],
                   "_pst_eff_a1_b0.5_phi",phi_ch[k],"_chr22.txt"))
prs = prs %>% select(V1,V2,V6)
prs.match = left_join(prs,summary.tar,by=c("V2"="rs_id"))
prs.file <- prs.match %>% 
  select(SNP,A1,BETA)
colSums(is.na(prs.file))
write.table(prs.file,file = paste0(temp.dir.prs,"prs_file"),col.names = T,row.names = F,quote=F)

res = system(paste0("/data/zhangh24/software/plink2_alpha --score-col-nums 3 --threads 2 --score ",temp.dir.prs,"prs_file  header no-mean-imputation --bfile ",temp.dir,"AFR_ref_chr22  --out ",temp.dir.prs,"prs_csx_",eth[i],"_phi",phi[k]))
time = proc.time()-time1
save(time,file = paste0(out.dir,"prscsx_trep_",t_rep,"_phi_",k,".rdata"))
#system(paste0("mv ",temp.dir.prs,"/prs_csx_",eth[i],"_rho_",l,"_size_",m,"_GA_",i1,"_phi",phi[v],".sscore ",out.dir.sum,eth[i],"/prscsx/"))

