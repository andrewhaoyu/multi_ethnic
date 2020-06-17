args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
k = as.numeric(args[[2]])
l = as.numeric(args[[3]])
#use plink2 to calculate prs
#load LD_clump_file

library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01)

n.snp.mat <- cbind(expand.grid(pthres,pthres),0,0,0)
colnames(n.snp.mat) = c("pthres1","pthres2","nsnp")
#load snp infor data to get maf to standarize effect size
load("/data/zhangh24/multi_ethnic/result/LD_simulation/snp.infor.rdata")
colnames(snp.infor)[1] <- "SNP"
snp.infor.select = snp.infor %>% 
  select(SNP,EUR,AFR,AMR,EAS)


#read summary level statistics from EUR
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

#for(i in 2:length(eth)){
  #read LD clumped SNPs
  LD <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/LD_clump_two_dim.clumped")))
  clump.snp <- LD[,3,drop=F] 
  #read the target ethnic group summary level statistics
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
  #k for the p-value threshold on the target population
  #l for the p-value threshold on the eur population
 # for(k in 1:length(pthres)){
   # for(l in 1:length(pthres)){
      print(c(i,k,l))
      #select the SNPs passing the threshold
      prs.select = prs.all %>% filter(P<=pthres[k]|
                           peur<=pthres[l]) 
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
      write.table(prs.file,file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_file_pvalue_eb_",k,"_",l),col.names = T,row.names = F,quote=F)  
    #}
    
  #}
  
#}
#write.csv(n.snp.mat,file =paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/n_snp_mat.csv"))

#use plink score funciton to get prs
code <- rep("c",10000)
temp <- 1
for(i in 2:length(eth)){
  for(j in 1:22){
    for(k in 1:length(pthres)){
      for(l in 1:length(pthres)){
        temp.code <- paste0("/data/zhangh24/software/plink2 --score /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/prs_file_pvalue_eb_",k,"_",l," no-sum no-mean-imputation  --allow-no-sex --bfile /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/prs/chr",j,"_prs_eb_",k,"_",l)
        code[temp] <- temp.code
        temp <- temp+1
      }
    }
  }

}
code <- code[1:(temp-1)]
write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/calculate_prs_eb.sh"),col.names = F,row.names = F,quote=F)
# 
# 
# write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/test.sh"),col.names = F,row.names = F,quote=F)




