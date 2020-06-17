#use plink2 to calculate for EUR best PRS, use EUR cofficients or AFR cofficients
#load LD_clump_file

library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
n.snp.mat <- matrix(0,length(pthres),4)
#load EUR r2 result
load("/data/zhangh24/multi_ethnic/result/LD_simulation/r2.mat.eur.rdata")
#get the best performance eur prs p-value threshold
jdx = which.max(r2.mat[,1])
i= 1
LD <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/LD_clump.clumped")))
clump.snp <- LD[,3,drop=F]  
sum.data <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out")))
colnames(sum.data)[2] <- "SNP"
prs.all <- left_join(clump.snp,sum.data,by="SNP") 
n.snp.mat <- matrix(0,1,4)
#use EUR regression coefficients
for(i in 2:length(eth)){
  
    prs.file <- prs.all %>% filter(P<=pthres[jdx]) %>% 
      select(SNP,A1,BETA)
    SNP = prs.file[,1,drop=F]
    save(SNP,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/EUR_best_PRS.rdata")
    #n.snp.mat[1,i] <- nrow(prs.file)
    write.table(prs.file,file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_eursnp_eurcoef"),col.names = T,row.names = F,quote=F)
  
  
}
#use target population regresssion coefficients
#get the number of SNPs in different population
n.snp.total = rep(0,4)
n.snp.total[1] = nrow(prs.file)
for(i in 2:length(eth)){
  prs.file <- prs.all %>% filter(P<=pthres[jdx]) %>% 
    select(SNP)
  sum.data <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out")))
  #combine the prs with the target population coefficients
  prs.file.com <- left_join(prs.file,sum.data,by="SNP")
  prs.file <- prs.file.com %>% select(SNP,A1,BETA)
  prs.file <- prs.file[complete.cases(prs.file),]
  n.snp.total[i] = nrow(prs.file)
  #n.snp.mat[j,i] <- nrow(prs.file)
  write.table(prs.file,file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_eursnp_tarcoef"),col.names = T,row.names = F,quote=F)
  
  
}


#use plink score funciton to get prs
code <- rep("c",10000)
temp <- 1
for(i in 2:length(eth)){
  for(j in 1:22){
    temp.code <- paste0("/data/zhangh24/software/plink2 --score /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_eursnp_eurcoef no-sum no-mean-imputation --bfile /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/chr",j,"_prs_eursnp_eurcoef")
    code[temp] <- temp.code
    temp <- temp+1
    
  }
  
}


for(i in 2:length(eth)){
  for(j in 1:22){
    temp.code <- paste0("/data/zhangh24/software/plink2 --score /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_eursnp_tarcoef no-sum no-mean-imputation --bfile /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/chr",j,"_prs_eursnp_tarcoef")
    code[temp] <- temp.code
    temp <- temp+1
    
  }
  
}

code <- code[1:(temp-1)]
write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/calculate_prs_eursnp.sh"),col.names = F,row.names = F,quote=F)



















#load snp infor data to get maf to standarize effect size
load("/data/zhangh24/multi_ethnic/result/LD_simulation/snp.infor.rdata")
colnames(snp.infor)[1] <- "SNP"
snp.infor.select = snp.infor %>% 
  select(SNP,EUR,AFR,AMR,EAS)

#use EB regresssion coefficients
#get the number of SNPs in different population
summary.eur <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[1],"/summary.out"),header=T))
colnames(summary.eur)[9] = "peur"
colnames(summary.eur)[7] = "beta_eur"
summary.eur = left_join(summary.eur,snp.infor.select,by="SNP")
summary.eur.select = summary.eur %>% 
  mutate(sd_eur = beta_eur/STAT,
         beta_eur_st=beta_eur*sqrt(2*EUR*(1-EUR)),
         sd_eur_st = sd_eur*sqrt(2*EUR*(1-EUR))) %>%
  rename(A1_eur = A1) %>% 
  select(SNP,beta_eur_st,sd_eur_st,peur,A1_eur)

load("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/EUR_best_PRS.rdata")
clump.snp <- SNP

#read the target ethnic group summary level statistics
for(i in 2:4){
  sum.data <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out")))
  colnames(sum.data)[2] <- "SNP"
  #get standarize the effect and sd
  sum.data.all = left_join(sum.data,snp.infor.select,by="SNP")
  
  sum.data.all = sum.data.all %>% 
    mutate(beta_st=BETA*sqrt(2*.[[9+i]]*(1-.[[9+i]])),
           beta_sd=BETA/STAT,
           sd_st=beta_sd*sqrt(2*.[[9+i]]*(1-.[[9+i]])))
  #combine the target level summary stat with EUR
  summary.com <- left_join(sum.data.all,summary.eur.select,by="SNP")
  #align EUR SNP effect to target population
  idx <- which(summary.com$A1!=summary.com$A1_eur)
  summary.com$beta_eur_st[idx] = -summary.com$beta_eur_st[idx]
  #combine the statistics with SNPs after clumping
  prs.all <- left_join(clump.snp,summary.com,by="SNP") 
  prs.select = prs.all
  #get the prior estimate using all the SNPs effect size
  prior.sigma = cov(cbind(prs.select$beta_st,
                          prs.select$beta_eur_st),use="complete.obs")
  load("/data/zhangh24/multi_ethnic/result/LD_simulation/causal_Sigma.rdata")
  true_sigma = Sigma
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
  post_beta_target = prs.select$BETA
  for(m in 1:nrow(prs.select)){
    Sigma = diag(c(prs.select$beta_sd[m]^2,
                   prs.select$sd_eur_st[m]^2))
    beta = c(prs.select$beta_st[m],
             prs.select$beta_eur_st[m])
    MAF.train = prs.select[m,9+i]
    MAF.ref = prs.select$EUR[m]
    
    if(is.na(Sigma[2,2])==0){
      post_beta_target[m] = PostBeta(beta,Sigma,prior.sigma,MAF.train,MAF.ref)[1] 
    }
  }
  prs.select$BETA_post = post_beta_target
  prs.file <- prs.select %>%
    select(SNP,A1,BETA_post)
  jdx = which(n.snp.mat[,1]==pthres[k]&n.snp.mat[,2]==pthres[l])
  # n.snp.mat[jdx,i+1] = nrow(prs.file)
  colnames(prs.file)[3] = "BETA"
  write.table(prs.file,file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_file_eursnp_eb"),col.names = T,row.names = F,quote=F)    
}


#use plink score funciton to get prs
code <- rep("c",10000)
temp <- 1
for(i in 2:length(eth)){
  for(j in 1:22){
    temp.code <- paste0("/data/zhangh24/software/plink2 --score /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_file_eursnp_eb no-sum no-mean-imputation --bfile /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/chr",j,"_prs_file_eursnp_eb")
    code[temp] <- temp.code
    temp <- temp+1
    
  }
  
}


code <- code[1:(temp-1)]
write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/calculate_prs_eursnp_eb.sh"),col.names = F,row.names = F,quote=F)




