#use plink2 to calculate prs
#load LD_clump_file
#i is the ethnic
#l is the causal proportion
#m is the training sample size
#k is the p value thres
#j is the number of chromsome
args = commandArgs(trailingOnly = T)
i_rep = as.numeric(args[[1]])
i = as.numeric(args[[2]])
j = as.numeric(args[[3]])
l = as.numeric(args[[4]])
m = as.numeric(args[[5]])
i1 = as.numeric(args[[6]])
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
#j = as.numeric(args[[3]])
sid <- Sys.getenv("SLURM_JOB_ID")
dir.create(paste0("/lscratch/",sid,"/",eth[i],"/"),showWarnings = F)
temp.dir <- paste0("/lscratch/",sid,"/",eth[i],"/")

system(paste0("cp ",cur.dir,eth[i],"/chr",j,".tag.* ",temp.dir))
system(paste0("ls ",temp.dir))

library(dplyr)
library(data.table)
setwd("/data/zhangh24/multi_ethnic/")

pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
#n.snp.mat <- matrix(0,length(pthres),4)
load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/snp.infor.rdata")
colnames(snp.infor)[1] <- "SNP"
snp.infor.select = snp.infor %>% 
  select(SNP,EUR,AFR,AMR,EAS,SAS)

#idx <- which(summary.eur$SNP=="rs2469601:80294998:A:G")
#for(l in 1:3){
  print(l)
 # for(m in 1:4){
    
    #read LD clumped SNPs
  LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".clumped")))
    clump.snp <- LD[,3,drop=F] 
    
    summary.data <- as.data.frame(fread(paste0(out.dir,eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))  
    colnames(sum.data)[2] <- "SNP"
    prs.all <- left_join(clump.snp,sum.data,by="SNP") 
    prs.file <- prs.all %>% filter(CHR==j) %>% 
      select(SNP,A1,BETA,STAT,P) %>% 
      rename(peur=P)
    summary.eur <- prs.file
    prs.file <- prs.all %>% filter(CHR==j) %>% 
      select(SNP)
    sum.data <- as.data.frame(fread(paste0(out.dir,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
    #combine the prs with the target population coefficients
    prs.file.com <- left_join(prs.file,sum.data,by="SNP")
    prs.file <- prs.file.com %>% select(SNP,A1,BETA,STAT,P)
    summary.tar <- prs.file
    
    idx <- which(summary.tar$A1!=summary.eur$A1)
    summary.eur$A1[idx] <- summary.tar$A1[idx]
    summary.eur$BETA[idx] <- -summary.eur$BETA[idx]
    
    summary.eur = left_join(summary.eur,snp.infor.select,by="SNP")
    summary.eur.select = summary.eur %>%
      mutate(sd_eur = BETA/STAT,
             beta_st=BETA*sqrt(2*EUR*(1-EUR)),
             sd_st = sd_eur*sqrt(2*EUR*(1-EUR))) %>%
      rename(MAF=EUR) %>% 
      select(SNP,beta_st,sd_st,A1,MAF,peur)
    
    summary.tar = left_join(summary.tar,snp.infor.select,by="SNP")
    summary.tar.select = summary.tar %>% 
      rename(MAF = eth[i]) %>% 
      mutate(sd_tar = BETA/STAT,
             beta_st=BETA*sqrt(2*MAF*(1-MAF)),
             sd_st = sd_tar*sqrt(2*MAF*(1-MAF))) %>%
      select(SNP,beta_st,sd_st,A1,MAF)
    #read the target ethnic group summary level statistics
    prior.sigma = cov(cbind(summary.tar.select$beta_st,
                            summary.eur.select$beta_st),use="complete.obs")
    #load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/causal_Sigma.rdata")
    #prior.sigma = Sigma
    #prior.sigma = (Sigma[c(i,1),c(i,1)])
    #implement the emprical Bayes
    #train is the target population
    #ref is EUR population
    PostBeta <- function(beta,Sigma,Sigma0,MAF.train,
                         MAF.ref){
      n <- length(beta)
      beta_post <- solve(solve(Sigma)+solve(Sigma0))%*%(solve(Sigma)%*%beta)
      beta_post[1] <- beta_post[1]/sqrt(2*MAF.train*(1-MAF.train))
      beta_post[2] <- beta_post[2]/sqrt(2*MAF.ref*(1-MAF.ref))
      return(beta_post)
    }
    post_beta_target = summary.tar.select$beta_st
    for(m_i in 1:nrow(summary.tar.select)){
      Sigma = diag(c(summary.tar.select$sd_st[m_i]^2,
                     summary.eur.select$sd_st[m_i]^2))
      beta = c(summary.tar.select$beta_st[m_i],
               summary.eur.select$beta_st[m_i])
      MAF.train = summary.tar.select$MAF[m_i]
      MAF.ref = summary.eur.select$MAF[m_i]
      
      if(is.na(det(Sigma))==0){
        if(det(Sigma)!=0){
          post_beta_target[m_i] = PostBeta(beta,Sigma,prior.sigma,MAF.train,MAF.ref)[1]   
        }
        
      }
    }
    summary.tar.select$BETA = post_beta_target
    summary.tar.select$P = summary.tar$P
    summary.tar.select$peur = summary.eur$peur
    prs.all = summary.tar.select
    for(k1 in 1:length(pthres)){
      for(k2 in 1:length(pthres)){
        
        prs.file <- prs.all %>% filter((P<=pthres[k1]|
                                          peur<=pthres[k2])) %>%
          select(SNP,A1,BETA)
        #setwd(temp.dir)
        
        #}
        #for(j in 1:22){
        if(nrow(prs.file)>0){
          write.table(prs.file,file = paste0(temp.dir,"prs_pvalue_eb_",k1,"_",k2,"_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1),col.names = T,row.names = F,quote=F)
          res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --score ",temp.dir,"prs_pvalue_eb_",k1,"_",k2,"_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1," no-sum no-mean-imputation --bfile ",temp.dir,"chr",j,".tag --exclude ",out.dir,eth[i],"/duplicated.id  --out ",out.dir,eth[i],"/prs/prs_eb_",k1,"_",k2,"_rho_",l,"_size_",m,"_chr_",j,"_rep_",i_rep,"_GA_",i1))
          if(res==2){
            stop()
          }
          #system(paste0("rm ",out.dir,eth[i],"prs/prs_eb_",l,"_size_",m,"_",j,"_rep_",i_rep,"_GA_",i1,".nosex"))
          #system(paste0("rm ",out.dir,"prs/prs_eb_",l,"_size_",m,"_",j,"_rep_",i_rep,"_GA_",i1,".log"))
          #system(paste0("rm ",out.dir,"prs/prs_eb_",l,"_size_",m,"_",j,"_rep_",i_rep,"_GA_",i1,".nopred"))
          system(paste0("rm ",temp.dir,"prs_pvalue_eb_",k1,"_",k2,"_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
          gc()
        }
        
      }
      
    }
    
    
  #}
#}

#pthres <- c(1E-10,1E-09,5E-08,1E-07,2.5E-07,5E-07,7.5E-07,1E-06,2.5E-06,5E-06,7.5E-06,1E-05,2.5e-05,5E-05,7.5e-05,1E-04,2.5E-04,5E-04,7.5E-04,1E-03)

system(paste0('rm -r /lscratch/',sid,'/',eth[i],'/'))

#}








