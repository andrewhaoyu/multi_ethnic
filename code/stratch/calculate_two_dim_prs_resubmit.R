#use plink2 to calculate prs
#load LD_clump_file
#i is the ethnic
#l is the causal proportion
#m is the training sample size
#k is the p value thres
#j is the number of chromsome
args = commandArgs(trailingOnly = T)
ind = as.numeric(args[[1]])
load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/mis_mat.rdata")
library(bc2,lib.loc = "/home/zhangh24/R/x86_64-pc-linux-gnu-library/4.0/")
size = 1000
start.end = startend(nrow(mis_mat),1000,ind)
start = start.end[1]
end = start.end[2]
row_ind =1 
ind.vec = mis_mat[row_ind,]
i_temp = ind.vec[2]

j_temp = ind.vec[4]

eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
#j = as.numeric(args[[3]])
sid <- Sys.getenv("SLURM_JOB_ID")



for(row_ind in start:end){
  ind.vec = mis_mat[row_ind,]
  i = ind.vec[2]
  i_rep = ind.vec[3]
  j = ind.vec[4]
  l = ind.vec[5]
  m = ind.vec[6]
  k1 = ind.vec[7]
  k2 = ind.vec[8]
  if(i!=i_temp|j!=j_temp){
    system(paste0('rm -r /lscratch/',sid,'/',eth[i_temp],'/'))
    i_temp = i
    j_temp = j
    dir.create(paste0('/lscratch/',sid,'/',eth[i_temp],"/"),showWarnings = F)
    temp.dir <- paste0('/lscratch/',sid,'/',eth[i_temp],"/")
    system(paste0("cp ",cur.dir,eth[i_temp],"/chr",j_temp,".tag.* ",temp.dir,"."))
    system(paste0("ls ",temp.dir))
    }
  #j = as.numeric(args[[3]])
  
  library(dplyr)
  library(data.table)
  setwd("/data/zhangh24/multi_ethnic/")
  
  pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
  #n.snp.mat <- matrix(0,length(pthres),4)
  
  
  #for(l in 1:3){
    print(l)
   # for(m in 1:4){
      
      print(m)
      summary.eur <- as.data.frame(fread(paste0("./result/LD_simulation_new/",eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep)))  
      colnames(summary.eur)[9] = "peur"
      colnames(summary.eur)[7] = "beta_eur"
      summary.eur.select = summary.eur %>% 
        select(SNP,beta_eur,peur)
      
      #read LD clumped SNPs
      LD <- as.data.frame(fread(paste0(cur.dir,eth[i],"/LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,".clumped")))
      clump.snp <- LD[,3,drop=F] 
      #read the target ethnic group summary level statistics
      sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_new/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep))) 
      colnames(sum.data)[2] <- "SNP"
      #combine the target level summary stat with EUR
      summary.com <- left_join(sum.data,summary.eur.select,by="SNP")
      #combine the statistics with SNPs after clumping
      prs.all <- left_join(clump.snp,summary.com,by="SNP") 
     # for(k1 in 1:length(pthres)){
      #  for(k2 in 1:length(pthres)){
          
          prs.file <- prs.all %>% filter((P<=pthres[k1]|
                                            peur<=pthres[k2])&
                                           CHR==j) %>% 
            select(SNP,A1,BETA)
          #setwd(temp.dir)
          
          #}
          #for(j in 1:22){
          if(nrow(prs.file)>0){
            write.table(prs.file,file = paste0(temp.dir,"prs_pvalue_two_dim_",k1,"_",k2,"_rho_",l,"_size_",m,"_rep_",i_rep),col.names = T,row.names = F,quote=F)
            system(paste0("/data/zhangh24/software/plink2 --threads 2 --score ",temp.dir,"prs_pvalue_two_dim_",k1,"_",k2,"_rho_",l,"_size_",m,"_rep_",i_rep," no-sum no-mean-imputation --bfile ",temp.dir,"chr",j,".tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out ",cur.dir,eth[i],"/prs/prs_two_dim_",k1,"_",k2,"_rho_",l,"_size_",m,"_chr_",j,"_rep_",i_rep))
            system(paste0("rm ",cur.dir,eth[i],"/prs/prs_two_dim_",k1,"_",k2,"_rho_",l,"_size_",m,"_chr_",j,"_rep_",i_rep,".nosex"))
            system(paste0("rm ",cur.dir,eth[i],"/prs/prs_two_dim_",k1,"_",k2,"_rho_",l,"_size_",m,"_chr_",j,"_rep_",i_rep,".nopred"))
            system(paste0("rm ",cur.dir,eth[i],"/prs/prs_two_dim_",k1,"_",k2,"_rho_",l,"_size_",m,"_chr_",j,"_rep_",i_rep,".log"))
            #system(paste0("/data/zhangh24/software/plink2 --score ",cur.dir,eth[i],"/prs/prs_file_pvalue_",k,"_rho_",l,"_size_",m,,"_rep_",i_rep," no-sum no-mean-imputation --bfile ",cur.dir,eth[i],"/all_chr.tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out ",cur.dir,eth[i],"/prs/prs_",k,"_rho_",l,"_size_",m))
          }
          
       # }
        
      #}
      
      
  #  }
 # }
  
  #pthres <- c(1E-10,1E-09,5E-08,1E-07,2.5E-07,5E-07,7.5E-07,1E-06,2.5E-06,5E-06,7.5E-06,1E-05,2.5e-05,5E-05,7.5e-05,1E-04,2.5E-04,5E-04,7.5E-04,1E-03)
  
  
  gc()
}


#}








