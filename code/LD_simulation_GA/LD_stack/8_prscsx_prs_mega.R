#use plink2 to calculate prs
#load LD_clump_file
#i is the ethnic
#l is the causal proportion
#m is the training sample size
#k is the p value thres
#j is the number of chromsome
#args = commandArgs(trailingOnly = T)

# i = as.numeric(args[[1]])
# l = as.numeric(args[[2]])
# m = as.numeric(args[[3]])
# i1 = as.numeric(args[[4]])
i = 2
l = 1
m = 1
#i_rep = as.numeric(args[[1]])
i1 = 1
library(dplyr)
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
eth <- c("EUR","AFR","AMR","EAS","SAS")

sid <- Sys.getenv("SLURM_JOB_ID")
dir.create(paste0('/lscratch/',sid,'/test/'),showWarnings = F)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/prs/'),showWarnings = FALSE)
temp.dir.prs = paste0('/lscratch/',sid,'/test/prs/')
library(data.table)
system(paste0("cp ",cur.dir,eth[i],"/all_chr_test.mega.* ",temp.dir))
system(paste0("ls ",temp.dir))
setwd("/data/zhangh24/multi_ethnic/")
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
setwd(paste0(out.dir.sum,eth[i],"/prscsx_mega/"))
phi = c("1","1e-02","1e-04","1e-06")
for(v in 1:4){
  i_rep = 1
  #load target population posterior
  prs = fread(paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",eth[i],
                     "_pst_eff_a1_b0.5_phi",phi[v],".txt"))
  prs_infor = prs %>% 
    select(V2,V4)
  colnames(prs_infor) = c("rs_id","A1")
  n_rep = 10
  Beta = matrix(0,nrow(prs),n_rep)
  for(i_rep in 1:n_rep){
    print(i_rep)
    prs = fread(paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",eth[i],
                       "_pst_eff_a1_b0.5_phi",phi[v],".txt"))
    Beta[,i_rep] = prs$V6
    
  }
  
  prs = cbind(prs_infor,Beta)
  #colnames(prs) = c("CHR","rs_id","POS","A1","A2","BETA")
  
  snp.infor = snp.infor.match %>% 
    rename(SNP=id) %>% 
    select(SNP,rs_id)
  
  prs.match = left_join(prs,snp.infor,by=c("rs_id"="rs_id"))
  
  prs.file <- prs.match %>% 
    select(SNP,A1,paste0("V",c(1:n_rep)))
  colSums(is.na(prs.file))
  write.table(prs.file,file = paste0(temp.dir.prs,"prs_file"),col.names = T,row.names = F,quote=F)
  
  res = system(paste0("/data/zhangh24/software/plink2_alpha --score-col-nums 3,4,5,6,7,8,9,10,11,12 --threads 2 --score ",temp.dir.prs,"prs_file  cols=+scoresums,-scoreavgs header no-mean-imputation --bfile ",temp.dir,"all_chr_test.mega  --out ",temp.dir.prs,"prs_csx_",eth[i],"_rho_",l,"_size_",m,"_GA_",i1,"_phi",phi[v]))
  
  system(paste0("mv ",temp.dir.prs,"/prs_csx_",eth[i],"_rho_",l,"_size_",m,"_GA_",i1,"_phi",phi[v],".sscore ",out.dir.sum,eth[i],"/prscsx_mega/"))
  
  
  #load eur population posterior
  i_rep = 1
  #load the first element for the snp information
  prs = fread(paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_EUR_pst_eff_a1_b0.5_phi"
                     ,phi[v],".txt"))
  prs_infor = prs %>% 
    select(V2,V4)
  colnames(prs_infor) = c("rs_id","A1")
  n_rep = 10
  Beta = matrix(0,nrow(prs),n_rep)
  for(i_rep in 1:n_rep){
    print(i_rep)
    prs = fread(paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_EUR_pst_eff_a1_b0.5_phi"
                       ,phi[v],".txt"))
    Beta[,i_rep] = prs$V6
    
  }
  
  prs = cbind(prs_infor,Beta)
  #colnames(prs) = c("CHR","rs_id","POS","A1","A2","BETA")
  
  snp.infor = snp.infor.match %>% 
    rename(SNP=id) %>% 
    select(SNP,rs_id)
  
  prs.match = left_join(prs,snp.infor,by=c("rs_id"="rs_id"))
  
  prs.file <- prs.match %>% 
    select(SNP,A1,paste0("V",c(1:n_rep)))
  colSums(is.na(prs.file))
  write.table(prs.file,file = paste0(temp.dir.prs,"prs_file"),col.names = T,row.names = F,quote=F)
  
  res = system(paste0("/data/zhangh24/software/plink2_alpha --score-col-nums 3,4,5,6,7,8,9,10,11,12 --threads 2 --score ",temp.dir.prs,"prs_file cols=+scoresums,-scoreavgs header no-mean-imputation --bfile ",temp.dir,"all_chr_test.mega  --out ",temp.dir.prs,"prs_csx_EUR_rho_",l,"_size_",m,"_GA_",i1,"_phi",phi[v]))
  system(paste0("mv ",temp.dir.prs,"/prs_csx_EUR_rho_",l,"_size_",m,"_GA_",i1,"_phi",phi[v],".sscore ",out.dir.sum,eth[i],"/prscsx_mega/"))
  
}





#temp = fread(paste0(out.dir.sum,eth[i],"/prscsx/prs_csx_rho_",l,"_size_",m,"_GA_",i1,".sscore"))
