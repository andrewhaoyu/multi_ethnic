#select SNPs based on best EUR prs, use coefficients of EUR population
#i represent ethnic group
#l represent causal proportion
#m represent sample size
args = commandArgs(trailingOnly = T)
i_rep = as.numeric(args[[1]])
i = as.numeric(args[[2]])
i1 = as.numeric(args[[3]])
setwd("/data/zhangh24/multi_ethnic/")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
#where the summary data located
out.dir.sum <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)

sid <- Sys.getenv("SLURM_JOB_ID")
#system(paste0("rm -rf /lscratch/",sid,"/",eth[i]))
dir.create(paste0("/lscratch/",sid,"/",eth[i],"/"),showWarnings = F)
temp.dir <- paste0("/lscratch/",sid,"/",eth[i],"/")
system(paste0("cp ",cur.dir,eth[i],"/all_chr_test.mega.* ",temp.dir,"."))
system(paste0("ls ",temp.dir))
library(dplyr)
library(data.table)

for(l in 1:3){
  for(m in 1:4){
    #n.snp.mat <- matrix(0,length(pthres),4)
    #load EUR r2 result
    load(paste0(out.dir,"LD.clump.result.CT.rdata"))
    
    #load("/data/zhangh24/multi_ethnic/result/LD_simulation/r2.mat.eur.rdata")
    #keep the EUR results with sample size at 100,000 (m = 4)
    r2.mat <- as.data.frame(LD.result.list[[2]]) %>% 
      filter(eth_vec=="EUR"&
               m_vec==4&
               l_vec==l&
               ga_vec==i1) %>% 
      filter(r2.vec.test==max(r2.vec.test))
    #get the best performance eur prs p-value threshold
    p.cut = r2.mat$pthres_vec
    r_ind = 1
    w_ind = 1
    LD <- as.data.frame(fread(paste0(out.dir,eth[1],"/LD_clump_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    
    clump.snp <- LD[,3,drop=F] 
    
    sum.data <- as.data.frame(fread(paste0(out.dir.sum,eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))  
    colnames(sum.data)[2] <- "SNP"
    prs.all <- left_join(clump.snp,sum.data,by="SNP") 
    prs.file <- prs.all %>% filter(P<=p.cut) %>% 
      select(SNP,A1,BETA)
    prs.eur = prs.file
    
    #load the target ethnic coefficients
    sum.tar <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
    sum.tar = sum.tar %>% 
      mutate(A1_tar = A1,
             BETA_tar = BETA) %>% 
      select(SNP,A1_tar,BETA_tar)
    prs.eur = prs.eur %>% 
      mutate(A1_eur = A1,
             BETA_eur = BETA)
    
    prs.file = left_join(prs.eur,sum.tar,by="SNP") 
   
    prs.file = prs.file %>% 
      mutate(BETA_tar = ifelse(A1_tar==A1,BETA_tar,-BETA_tar),
        A1_tar = ifelse(A1_tar==A1,A1_tar,A1_eur))
    all.equal(prs.file[,"A1_tar"],
              prs.file[,"A1_eur"])
    
    
    #align the two different ethnic groups
    idx <- which(prs.tar$A1!=prs.eur$A1)
    if(length(idx)>0){
      prs.eur$A1[idx] = prs.tar$A1[idx]
      prs.eur$BETA[idx] = -prs.eur$BETA[idx]
    }
    
    
    
    
    
    
    write.table(prs.file,file = paste0("/lscratch/",sid,"/",eth[i],"/prs_eursnp_eurcoef_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1),col.names = T,row.names = F,quote=F)
    res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --score /lscratch/",sid,'/',eth[i],"/prs_eursnp_eurcoef_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1, " no-sum no-mean-imputation --bfile ",temp.dir,"/chr",j,".tag --exclude ",out.dir,eth[i],"/duplicated.id  --out ",out.dir,eth[i],"/prs/prs_eursnp_eurcoef_rho_",l,"_size_",m,"_",j,"_rep_",i_rep,"_GA_",i1))
    if(res==2){
      stop()
    }
    system(paste0("rm ",out.dir,eth[i],"/prs/prs_eursnp_eurcoef_rho_",l,"_size_",m,"_",j,"_rep_",i_rep,"_GA_",i1,".nosex"))
    system(paste0("rm ",out.dir,eth[i],"/prs/prs_eursnp_eurcoef_rho_",l,"_size_",m,"_",j,"_rep_",i_rep,"_GA_",i1,".log"))
    system(paste0("rm ",out.dir,eth[i],"/prs/prs_eursnp_eurcoef_rho_",l,"_size_",m,"_",j,"_rep_",i_rep,"_GA_",i1,".nopred"))
    system(paste0("rm /lscratch/",sid,"/",eth[i],"/prs_eursnp_eurcoef_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
    gc()
    #system(paste0('ls/lscratch/',sid,"/"))
  }
}
system(paste0('rm -rf /lscratch/',sid,'/',eth[i],'/'))
