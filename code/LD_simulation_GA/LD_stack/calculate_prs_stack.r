#use plink2 to calculate prs
#load LD_clump_file
#i is the ethnic
#l is the causal proportion
#m is the training sample size
#k is the p value thres
#j is the number of chromsome
args = commandArgs(trailingOnly = T)

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
i_rep = as.numeric(args[[3]])
j = as.numeric(args[[4]])
#l = as.numeric(args[[3]])
#m = as.numeric(args[[4]])
i1 = as.numeric(args[[5]])
#i_rep = 1

#j = as.numeric(args[[3]])
library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
#n.snp.mat <- matrix(0,length(pthres),4)
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
sid <- Sys.getenv("SLURM_JOB_ID")
dir.create(paste0('/lscratch/',sid,'/test/'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/prs/'),showWarnings = FALSE)
temp.dir.prs = paste0('/lscratch/',sid,'/test/prs/')
system(paste0("cp ",cur.dir,eth[i],"/chr",j,".mega.* ",temp.dir,"."))
system(paste0("ls ",temp.dir))
print("step1 finished")
#for(l in 1:3){
r2_vec = c(0.01,0.05,0.1,0.2,0.5)
wc_base_vec = c(50,100,200,500)
setwd("/data/zhangh24/multi_ethnic/")  
for(m in 1:4){
  sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
  for(r_ind in 1:length(r2_vec)){
    wc_vec = wc_base_vec/r2_vec[r_ind]
    for(w_ind in 1:length(wc_vec)){
      print(c(r_ind,w_ind))
      
  
  
  LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
  clump.snp <- LD[,3,drop=F] 
  #here use summary_out or summary_MAF_out doesn't make difference
  #since the LD_clumping step only include MAF_out
  
  
  colnames(sum.data)[2] <- "SNP"
  prs.all <- left_join(clump.snp,sum.data,by="SNP")
  n_pthres <- length(pthres)
  q_range = data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres))
  
  prs.file <- prs.all %>% filter(CHR==j) %>% 
    select(SNP,A1,BETA)
  write.table(prs.file,file = paste0(temp.dir.prs,"prs_file"),col.names = T,row.names = F,quote=F)
  p.value.file <- prs.all %>% filter(CHR==j) %>% 
    select(SNP,P)
  write.table(p.value.file,file = paste0(temp.dir.prs,"p_value_file"),col.names = T,row.names = F,quote=F)
  
  temp = 1
  for(k in 1:length(pthres)){
    prs.file <- prs.all %>% filter(P<=pthres[k]&CHR==j) 
    #  select(SNP,A1,BETA)
  
    
    if(nrow(prs.file)>0){
      q_range[temp,1] = paste0("p_value_",k)
      q_range[temp,3] = pthres[k]
      temp = temp+1
     }
    
    
  }
  q_range = q_range[1:(temp-1),]
  write.table(q_range,file = paste0(temp.dir.prs,"q_range_file"),row.names = F,col.names = F,quote=F)
  res <- system(paste0("/data/zhangh24/software/plink2 --q-score-range ",temp.dir.prs,"q_range_file ",temp.dir.prs,"p_value_file header --threads 2 --score ",temp.dir.prs,"prs_file header no-sum no-mean-imputation --bfile ",temp.dir,"chr",j,".mega --exclude /data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/duplicated.id  --out ",temp.dir.prs,"prs_chr_",j,"_rho_",l,"_size_",m,"_chr_",j,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind))
  #system(paste0("ls ",temp.dir.prs))
  #system(paste0("/data/zhangh24/software/plink2 --score ",cur.dir,eth[i],"/prs/prs_file_pvalue_",k,"_rho_",l,"_size_",m,,"_rep_",i_rep," no-sum no-mean-imputation --bfile ",cur.dir,eth[i],"/all_chr.tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out ",cur.dir,eth[i],"/prs/prs_",k,"_rho_",l,"_size_",m))
  print("step2 finished")
  if(res==2){
    stop()
  }
  
  res = system(paste0("mv ",temp.dir.prs,"*.profile ",out.dir,eth[i],"/prs/"))
  if(res==2){
    stop()
  }
  system(paste0("rm -rf ",temp.dir.prs))
  dir.create(paste0(temp.dir.prs),showWarnings = FALSE)
  
}
  }
}
#}
print("step3 finished")
system(paste0('rm -r /lscratch/',sid,'/',eth[i],'/'))







# result.matrix <- matrix(0,3,1)
# for(l in 1:3){
#   #for(i1 in 1:2){
#     sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_new/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep)))
#     idx <- which(sum.data$P<=5E-08)
# 
#     result.matrix[l,i1] <- length(idx)
#   #}
# }
