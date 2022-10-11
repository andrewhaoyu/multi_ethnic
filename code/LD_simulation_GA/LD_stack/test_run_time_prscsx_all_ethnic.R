args = commandArgs(trailingOnly = T)
t_rep = as.numeric(args[[1]])
k = as.numeric(args[[2]])
#test the running time and memory for TDLD-SLEB and PRS-CSx on CHR 22
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/1KGLD'),showWarnings = FALSE)

system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_eur ",temp.dir,"1KGLD"))
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_afr ",temp.dir,"1KGLD"))
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_amr ",temp.dir,"1KGLD"))
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_eas ",temp.dir,"1KGLD"))
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_sas ",temp.dir,"1KGLD"))
system(paste0("cp /data/zhangh24/software/PRScsx/1KGLD/snpinfo_mult_1kg_hm3 ",temp.dir,"1KGLD"))

time1 = proc.time()
i = 2
phi = c(1,1E-02,1E-04,1E-06)
library(data.table)
library(tidyverse)
eth <- c("EUR","AFR","AMR","EAS","SAS")
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
              " --sst_file=",path_to_sum,"EUR_sumstats.txt,",path_to_sum,"AFR_sumstats.txt,",path_to_sum,"AMR_sumstats.txt,",path_to_sum,"EAS_sumstats.txt,",path_to_sum,"SAS_sumstats.txt",
              " --n_gwas=100000,15000,15000,15000,15000",
              " --pop=EUR,AFR,AMR,EAS,SAS",
              " --chrom=",j,
              " --phi=",phi[k], 
              " --out_dir=",out.dir,
              " --out_name=test_five"))
#calculate prs for the two coefficients
phi_ch = c("1e+00","1e-02","1e-04","1e-06")
#calcualte prs for all
prs = fread(paste0(out.dir,"test_five_",eth[i],
                   "_pst_eff_a1_b0.5_phi",phi_ch[k],"_chr22.txt"))
prs = prs %>% select(V1,V2,V6)
#load AFR summary data
load(paste0(out.dir,"tar_sumdata.rdata"))
prs.match = left_join(prs,summary.tar,by=c("V2"="rs_id"))
dir.create(paste0('/lscratch/',sid,'/test/prs/'),showWarnings = FALSE)
temp.dir.prs = paste0('/lscratch/',sid,'/test/prs/')
system(paste0("cp ",out.dir,"AFR_test_mega_chr22.bed ",temp.dir,"AFR_test_mega_chr22.bed"))
system(paste0("cp ",out.dir,"AFR_test_mega_chr22.bim ",temp.dir,"AFR_test_mega_chr22.bim"))
system(paste0("cp ",out.dir,"AFR_test_mega_chr22.fam ",temp.dir,"AFR_test_mega_chr22.fam"))
time2 = proc.time()
data_list = list()
for(i_eth in 1:5){
  data = fread(paste0(out.dir,"test_five_",eth[i_eth],
                      "_pst_eff_a1_b0.5_phi",phi_ch[k],"_chr22.txt"))
  data_temp = data %>% select(V2, V4)
  data_list[[i_eth]] = data_temp
}
snp_set = rbindlist(data_list) %>% distinct()
colnames(snp_set) = c("SNP", "A1")
BETA_mat = matrix(0,nrow(snp_set),2)
temp = 1
for(i_eth in 1:5){
  data = fread(paste0(out.dir,"test_five_",eth[i_eth],
                      "_pst_eff_a1_b0.5_phi",phi_ch[k],"_chr22.txt"))
  data_tar = data %>% select(V2, V6)
  snp_set_temp = left_join(snp_set, data_tar, by = c("SNP"="V2")) %>% 
    mutate(BETA = ifelse(is.na(V6),0,V6))
  BETA_mat[, temp] = snp_set_temp$BETA
  temp = temp + 1
}
load(paste0(out.dir,"tar_sumdata.rdata"))
summary.tar.select = summary.tar %>% 
  select(SNP,rs_id)
prs_file = cbind(snp_set,BETA_mat)
colnames(prs_file)[1] = "rs_id"
prs.match = left_join(prs_file,summary.tar.select,by="rs_id")
#replace the rs_id with SNP ID in 1KG genomes
prs_file[,1] = prs.match[,5]
colnames(prs_file)[1] = "SNP"

n_col = ncol(prs_file)
#file_out = paste0("/data/zhangh24/multi_ethnic/result/AOU/prs/PRSCSx/",eth[i],"/",trait,"")
write.table(prs_file,file = paste0(temp.dir,"prs_file"),col.names = T,row.names = F,quote=F)

ncol = ncol(prs_file)
res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                    "--score-col-nums 3-",n_col," --threads 2 ",
                    "--score ",temp.dir,"prs_file cols=+scoresums,-scoreavgs header no-mean-imputation ",
                    "--bfile ", temp.dir,"AFR_test_mega_chr22 ",
                    "--out ",temp.dir.prs,"prs_csx_",eth[i],"_phi",phi[k]))

time = proc.time()-time1
time_prs = proc.time()-time2
time_list = c(time, time_prs)
save(time_list,file = paste0(out.dir,"prscsx_five_trep_",t_rep,"_phi_",k,".rdata"))
#system(paste0("mv ",temp.dir.prs,"/prs_csx_",eth[i],"_rho_",l,"_size_",m,"_GA_",i1,"_phi",phi[v],".sscore ",out.dir.sum,eth[i],"/prscsx/"))

