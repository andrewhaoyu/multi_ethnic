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
      mutate(SE_eur = BETA/STAT) %>% 
      select(SNP,A1,BETA,SE_eur)
    prs.eur = prs.file
    
    #load the target ethnic coefficients
    sum.tar <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
    sum.tar = sum.tar %>% 
      mutate(A1_tar = A1,
             BETA_tar = BETA,
             SE_tar = BETA/STAT) %>% 
      select(SNP,A1_tar,BETA_tar,SE_tar)
    prs.eur = prs.eur %>% 
      mutate(A1_eur = A1,
             BETA_eur = BETA)
    
    prs.file = left_join(prs.eur,sum.tar,by="SNP") 
   
    prs.file = prs.file %>% 
      mutate(BETA_tar = ifelse(A1_tar==A1,BETA_tar,-BETA_tar),
        A1_tar = ifelse(A1_tar==A1,A1_tar,A1_eur))
    all.equal(prs.file[,"A1_tar"],
              prs.file[,"A1_eur"])
    
    #load EB coefficients
    prs.file = prs.file %>% 
      mutate(z_stat_tar = BETA_tar/SE_tar,
             z_stat_eur = BETA_eur/SE_eur)
    z.mat = prs.file %>% 
      select(z_stat_tar,z_stat_eur)
    prior.sigma = cov(z.mat,use="complete.obs")-diag(2)
    post.sigma = solve(solve(prior.sigma)+diag(2))
    z.post = as.matrix(z.mat)%*%post.sigma
    post_beta_tar = z.post[,"z_stat_tar"]*prs.file$SE_tar
    
    prs.file = prs.file %>% 
      mutate(post_beta_tar=ifelse(is.na(post_beta_tar),BETA_tar,post_beta_tar))
    
    prs.coef = prs.file %>% 
      select(SNP,A1,BETA_eur,BETA_tar,post_beta_tar)
    
    write.table(prs.coef,file = paste0("/lscratch/",sid,"/",eth[i],"/prs_besteur_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1),col.names = T,row.names = F,quote=F)
    res = system(paste0("/data/zhangh24/software/plink2_alpha --score-col-nums 3,4,5 --threads 2 --score /lscratch/",sid,'/',eth[i],"/prs_besteur_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1, " header no-mean-imputation --bfile ",temp.dir,"/all_chr_test.mega --exclude ",out.dir.sum,eth[i],"/duplicated.id  --out /lscratch/",sid,'/',eth[i],"/prs_besteur_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
    system(paste0("mv /lscratch/",sid,'/',eth[i],"/prs_besteur_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".sscore ",out.dir,eth[i],"/prs/"))
    if(res==2){
      stop()
    }
    
    gc()
    #system(paste0('ls/lscratch/',sid,"/"))
  }
}
system(paste0('rm -rf /lscratch/',sid,'/',eth[i],'/'))
