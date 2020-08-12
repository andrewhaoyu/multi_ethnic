#select SNPs based on best EUR prs, use EB coefficients 
#i represent ethnic group
#l represent causal proportion
#m represent sample size
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
j = as.numeric(args[[4]])
library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01)
#n.snp.mat <- matrix(0,length(pthres),4)
#load EUR r2 result
setwd("/data/zhangh24/multi_ethnic/")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
load(paste0(cur.dir,"LD.clump.result.rdata"))
load("/data/zhangh24/multi_ethnic/result/LD_simulation/r2.mat.eur.rdata")
#keep the EUR results with sample size at 100,000 (m = 4)
r2.mat <- as.data.frame(LD.result.list[[2]]) %>% 
  filter(eth.vec=="EUR"&
           m_vec==4&
           l_vec==l)

#keep the EUR results with sample size at 100,000 (m = 4)
#get the best performance eur prs p-value threshold
k = which.max(r2.mat$r2.vec)
LD <- as.data.frame(fread(paste0(cur.dir,eth[1],"/LD_clump_rho_",l,"_size_",4,".clumped")))
clump.snp <- LD[,3,drop=F]  
sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_new/",eth[1],"/summary_out_rho_",l,"_size_",4)))  
colnames(sum.data)[2] <- "SNP"
prs.all <- left_join(clump.snp,sum.data,by="SNP") 
prs.file <- prs.all %>% filter(P<=pthres[k]) %>% 
  select(SNP,A1,BETA,STAT)
summary.eur <- prs.file
#use target population regresssion coefficients
#get the number of SNPs in different population
prs.file <- prs.all %>% filter(P<=pthres[k]) %>% 
  select(SNP)
sum.data <- as.data.frame(fread(paste0(cur.dir,eth[i],"/summary_out_rho_",l,"_size_",m)))
#combine the prs with the target population coefficients
prs.file.com <- left_join(prs.file,sum.data,by="SNP")
prs.file <- prs.file.com %>% select(SNP,A1,BETA,STAT)
summary.tar <- prs.file
#align the two population coefficients
idx <- which(summary.eur$A1!=summary.tar$A1)
summary.tar$A1[idx] <- summary.eur$A1[idx]
summary.tar$BETA[idx] <- -summary.tar$BETA[idx]

#load snp infor data to get maf to standarize effect size
load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/snp.infor.rdata")
colnames(snp.infor)[1] <- "SNP"
snp.infor.select = snp.infor %>% 
  select(SNP,EUR,AFR,AMR,EAS,SAS)

#use EB regresssion coefficients
#get the number of SNPs in different population
summary.eur = left_join(summary.eur,snp.infor.select,by="SNP")
summary.eur.select = summary.eur %>%
  mutate(sd_eur = BETA/STAT,
         beta_st=BETA*sqrt(2*EUR*(1-EUR)),
         sd_st = sd_eur*sqrt(2*EUR*(1-EUR))) %>%
  rename(MAF=EUR) %>% 
  select(SNP,beta_st,sd_st,A1,MAF)
summary.tar = left_join(summary.tar,snp.infor.select,by="SNP")
summary.tar.select = summary.tar %>% 
  rename(MAF = eth[i]) %>% 
  mutate(sd_tar = BETA/STAT,
         beta_st=BETA*sqrt(2*MAF*(1-MAF)),
         sd_st = sd_tar*sqrt(2*MAF*(1-MAF))) %>%
  select(SNP,beta_st,sd_st,A1,MAF)
#read the target ethnic group summary level statistics
  prior.sigma = cov(cbind(summary.eur.select$beta_st,
                          summary.tar.select$beta_st),use="complete.obs")
  load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/causal_Sigma.rdata")
  prior.sigma = Sigma
  prior.sigma = (Sigma[c(i,1),c(i,1)])
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
    
    if(is.na(Sigma[2,2])==0){
      post_beta_target[m_i] = PostBeta(beta,Sigma,prior.sigma,MAF.train,MAF.ref)[1] 
    }
  }
  summary.tar.select$BETA = post_beta_target
  
  prs.file <- summary.tar.select %>%
    select(SNP,A1,BETA)
  prs.file = prs.file[complete.cases(prs.file),]
  write.table(prs.file,file = paste0(cur.dir,eth[i],"/prs/prs_file_eursnp_eb_rho_",l,"_size_",m),col.names = T,row.names = F,quote=F)    

  system(paste0("/data/zhangh24/software/plink2 --threads 2 --score ",cur.dir,eth[i],"/prs/prs_file_eursnp_eb_rho_",l,"_size_",m," no-sum no-mean-imputation --bfile ",cur.dir,eth[i],"/chr",j,".tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out ",cur.dir,eth[i],"/prs/prs_eursnp_eb_rho_",l,"_size_",m,"_",j))




