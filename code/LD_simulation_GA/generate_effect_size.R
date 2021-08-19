#two different genetic architecture is used here
#the first one assumes same total heritability for all common causal SNPs
#the second one assumes same total heritability for all  causal SNPs
#generate effect size for the causal SNPs
#genetic correlation 0.8,0.6
#strong mild no negative selection
#two Genetic architecture, fixed heritbaility for common snps
#fixed heritibility for all snps
GenSigma <- function(sigma,n1,n2,n3,n4,n5,
                     gr12,gr13,gr14,gr15,
                     gr23,gr24,gr25,
                     gr34,gr35,
                     gr45){
  vecsigma <- sigma*c(1/n1,
                      gr12/sqrt(n1*n2),
                      gr13/sqrt(n1*n3),
                      gr14/sqrt(n1*n4),
                      gr15/sqrt(n1*n5),
                      gr12/sqrt(n1*n2),
                      1/n2,
                      gr23/sqrt(n2*n3),
                      gr24/sqrt(n2*n4),
                      gr25/sqrt(n2*n5),
                      gr13/sqrt(n1*n3),
                      gr23/sqrt(n2*n3),
                      1/n3,
                      gr34/sqrt(n3*n4),
                      gr35/sqrt(n3*n5),
                      gr14/sqrt(n1*n4),
                      gr24/sqrt(n2*n4),
                      gr34/sqrt(n3*n4),
                      1/n4,
                      gr45/sqrt(n4*n5),
                      gr15/sqrt(n1*n5),
                      gr25/sqrt(n2*n5),
                      gr35/sqrt(n3*n5),
                      gr45/sqrt(n4*n5),
                      1/n5)
  Sigma <- matrix(vecsigma,5,5)
  return(Sigma)
  
}

gr_vec = c(0.8,0.6)
alpha_vec = c(-1,-0.25,0)
genetic_arc_vec = c(1,2) 
load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/cau.snp.infor.list.rdata")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
#i for genetic correlation vec
#j for negative selection alphs
temp = 1
i1 = 1
# for(i in 1:2){
#   #for(j in 1:3){
#     for(l in 1:3){
      eth <- c("EUR","AFR","AMR","EAS","SAS")
      cau.snp.infor <- cau.snp.infor.list[[l]]
      library(data.table)
      #generate the effect-size for all the SNPs
      #even if the SNPs have 0 MAF in a particular population the effectsize will aslo be generated
      #it's okay since only the extracted genotype are those only existing in a particular population
      
      #get the number of causal SNPs in EUR
      EUR.idx <- which(cau.snp.infor$EUR>=0.01&
                         cau.snp.infor$EUR<=0.99)
      EUR.snp <- cau.snp.infor[EUR.idx,]
      n.EUR.snp <- nrow(EUR.snp)
      #get the number of causal SNPs in AFR
      AFR.idx <- which(cau.snp.infor$AFR>=0.01&
                         cau.snp.infor$AFR<=0.99)
      AFR.snp <- cau.snp.infor[AFR.idx,]
      n.AFR.snp <- nrow(AFR.snp)
      #get the number of causal SNPs in AMR
      AMR.idx <- which(cau.snp.infor$AMR>=0.01&
                         cau.snp.infor$AMR<=0.99)
      AMR.snp <- cau.snp.infor[AMR.idx,]
      n.AMR.snp <- nrow(AMR.snp)
      #get the number of causal SNPs in AMR
      EAS.idx <- which(cau.snp.infor$EAS>=0.01&
                         cau.snp.infor$EAS<=0.99)
      EAS.snp <- cau.snp.infor[EAS.idx,]
      n.EAS.snp <- nrow(EAS.snp)
      SAS.idx <- which(cau.snp.infor$SAS>=0.01&
                         cau.snp.infor$SAS<=0.99)
      SAS.snp <- cau.snp.infor[SAS.idx,]
      n.SAS.snp <- nrow(SAS.snp)
      
      n.total.snp <- nrow(cau.snp.infor)
      
      
      #first we need to generate the effect on standardized scale
      #
      sigma = 0.4
      gr12 = 0.8
      gr13 = 0.8
      gr14 = 0.8
      gr15 = 0.8
      gr23 = 0.8
      gr24 = 0.8
      gr25 = 0.8
      gr34 = 0.8
      gr35 = 0.8
      gr45 = 0.8
      n1 = as.numeric(n.EUR.snp)
      n2 = as.numeric(n.AFR.snp)
      n3 = as.numeric(n.AMR.snp)
      n4 = as.numeric(n.EAS.snp)
      n5 = as.numeric(n.SAS.snp)
      Sigma <- GenSigma(sigma,n1,n2,n3,n4,n5,
                        gr12,gr13,gr14,gr15,
                        gr23,gr24,gr25,
                        gr34,gr35,
                        gr45)
     # save(Sigma,file = paste0(cur.dir,"/causal_Sigma_",i1,".rdata"))
      library(mvtnorm)
      #beta represent standarize scale effect-size
      set.seed(666)
      beta =  rmvnorm(n.total.snp,c(0,0,0,0,0),
                      sigma=Sigma)
      colnames(beta) <- paste0("beta_",c("EUR","AFR","AMR","EAS","SAS"))
      #since this is the effect-size under standardized scale
      #to use GCTA, we need the original scale effect size
      EUR.bi <- cau.snp.infor$EUR>=0.01&
        cau.snp.infor$EUR<=0.99
      AFR.bi <-  cau.snp.infor$AFR>=0.01&
        cau.snp.infor$AFR<=0.99
      AMR.bi <-  cau.snp.infor$AMR>=0.01&
        cau.snp.infor$AMR<=0.99
      EAS.bi <-  cau.snp.infor$EAS>=0.01&
        cau.snp.infor$EAS<=0.99
      SAS.bi <-  cau.snp.infor$SAS>=0.01&
        cau.snp.infor$SAS<=0.99
      eth.bi <- cbind(EUR.bi,AFR.bi,AMR.bi,EAS.bi,SAS.bi)
      library(dplyr)
      MAF <- cau.snp.infor %>% select(EUR,AFR,AMR,EAS,SAS)
      MAF <- as.data.frame(MAF)
      for(i in 1:5){
        idx <- which(eth.bi[,i]==1)
        #beta.ori <- beta[idx,i]/sqrt(2*MAF[idx,i]*(1-MAF[idx,i]))
        select.cau <- cbind(cau.snp.infor[idx,1],
                            beta[idx,i])
        write.table(select.cau,file = paste0(cur.dir,eth[i],"/select.cau_rho",l,"_",i1),row.names = F,col.names = F,quote=F)
      }
      
      
#   }
# #}
# }
i1 = 1
#generate number of causal SNPs for l = 2
l = 2
n.cau = rep(0,5)
for(i in 1:5){
  select.cau <- read.table(paste0(cur.dir,eth[i],"/select.cau_rho",l,"_",i1),header=F)
  n.cau[i] = nrow(select.cau)
}
n.cau[c(2,3,4,1,5)]
select.cau1 = select.cau
herit <- nrow(select.cau)*var(select.cau$V2)
print(herit)


















#all SNPs causal, genetic correlation 0.8
i1 = 2
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
for(l in 1:3){
  eth <- c("EUR","AFR","AMR","EAS","SAS")
  cau.snp.infor <- cau.snp.infor.list[[l]]
  library(data.table)
  #generate the effect-size for all the SNPs
  #even if the SNPs have 0 MAF in a particular population the effectsize will aslo be generated
  #it's okay since only the extracted genotype are those only existing in a particular population
  
  #get the number of causal SNPs in EUR
  
  n.total.snp <- nrow(cau.snp.infor)
  
  
  sigma = 0.4
  gr12 = 0.8
  gr13 = 0.8
  gr14 = 0.8
  gr15 = 0.8
  gr23 = 0.8
  gr24 = 0.8
  gr25 = 0.8
  gr34 = 0.8
  gr35 = 0.8
  gr45 = 0.8
  n = as.numeric(n.total.snp)
  n1 = n
  n2 = n
  n3 = n
  n4 = n
  n5 = n
  Sigma <- GenSigma(sigma,n1,n2,n3,n4,n5,
                    gr12,gr13,gr14,gr15,
                    gr23,gr24,gr25,
                    gr34,gr35,
                    gr45)
  save(Sigma,file = paste0(cur.dir,"/causal_Sigma_",i1,".rdata"))
  library(mvtnorm)
  #beta represent standarize scale effect-size
  set.seed(666)
  beta =  rmvnorm(n.total.snp,c(0,0,0,0,0),
                  sigma=Sigma)
  colnames(beta) <- paste0("beta_",c("EUR","AFR","AMR","EAS","SAS"))
  #since this is the effect-size under standardized scale
  #to use GCTA, we need the original scale effect size
  EUR.bi <- cau.snp.infor$EUR>=0.01&
    cau.snp.infor$EUR<=0.99
  AFR.bi <-  cau.snp.infor$AFR>=0.01&
    cau.snp.infor$AFR<=0.99
  AMR.bi <-  cau.snp.infor$AMR>=0.01&
    cau.snp.infor$AMR<=0.99
  EAS.bi <-  cau.snp.infor$EAS>=0.01&
    cau.snp.infor$EAS<=0.99
  SAS.bi <-  cau.snp.infor$SAS>=0.01&
    cau.snp.infor$SAS<=0.99
  eth.bi <- cbind(EUR.bi,AFR.bi,AMR.bi,EAS.bi,SAS.bi)
  library(dplyr)
  MAF <- cau.snp.infor %>% select(EUR,AFR,AMR,EAS,SAS)
  MAF <- as.data.frame(MAF)
  total.herit = rep(0,5)
  for(i in 1:5){
    idx <- which(eth.bi[,i]==1)
    #beta.ori <- beta[idx,i]/sqrt(2*MAF[idx,i]*(1-MAF[idx,i]))
    select.cau <- cbind(cau.snp.infor[idx,1],
                        beta[idx,i])
  #  total.herit[i] = sum(select.cau$V2^2)
    write.table(select.cau,file = paste0(cur.dir,eth[i],"/select.cau_rho",l,"_",i1),row.names = F,col.names = F,quote=F)
  }
  
  total.herit[c(2,3,4,1,5)]
  
}









#genetic correlation 0.6, all SNPs causal
i1 = 3
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
for(l in 1:3){
  eth <- c("EUR","AFR","AMR","EAS","SAS")
  cau.snp.infor <- cau.snp.infor.list[[l]]
  library(data.table)
  #generate the effect-size for all the SNPs
  #even if the SNPs have 0 MAF in a particular population the effectsize will aslo be generated
  #it's okay since only the extracted genotype are those only existing in a particular population
  
  #get the number of causal SNPs in EUR
  
  n.total.snp <- nrow(cau.snp.infor)
  
  
  sigma = 0.4
  gr12 = 0.6
  gr13 = 0.6
  gr14 = 0.6
  gr15 = 0.6
  gr23 = 0.6
  gr24 = 0.6
  gr25 = 0.6
  gr34 = 0.6
  gr35 = 0.6
  gr45 = 0.6
  n = as.numeric(n.total.snp)
  n1 = n
  n2 = n
  n3 = n
  n4 = n
  n5 = n
  Sigma <- GenSigma(sigma,n1,n2,n3,n4,n5,
                    gr12,gr13,gr14,gr15,
                    gr23,gr24,gr25,
                    gr34,gr35,
                    gr45)
  save(Sigma,file = paste0(cur.dir,"/causal_Sigma_",i1,".rdata"))
  library(mvtnorm)
  #beta represent standarize scale effect-size
  set.seed(666)
  beta =  rmvnorm(n.total.snp,c(0,0,0,0,0),
                  sigma=Sigma)
  colnames(beta) <- paste0("beta_",c("EUR","AFR","AMR","EAS","SAS"))
  #since this is the effect-size under standardized scale
  #to use GCTA, we need the original scale effect size
  EUR.bi <- cau.snp.infor$EUR>=0.01&
    cau.snp.infor$EUR<=0.99
  AFR.bi <-  cau.snp.infor$AFR>=0.01&
    cau.snp.infor$AFR<=0.99
  AMR.bi <-  cau.snp.infor$AMR>=0.01&
    cau.snp.infor$AMR<=0.99
  EAS.bi <-  cau.snp.infor$EAS>=0.01&
    cau.snp.infor$EAS<=0.99
  SAS.bi <-  cau.snp.infor$SAS>=0.01&
    cau.snp.infor$SAS<=0.99
  eth.bi <- cbind(EUR.bi,AFR.bi,AMR.bi,EAS.bi,SAS.bi)
  library(dplyr)
  MAF <- cau.snp.infor %>% select(EUR,AFR,AMR,EAS,SAS)
  MAF <- as.data.frame(MAF)
  for(i in 1:5){
    idx <- which(eth.bi[,i]==1)
    #beta.ori <- beta[idx,i]/sqrt(2*MAF[idx,i]*(1-MAF[idx,i]))
    select.cau <- cbind(cau.snp.infor[idx,1],
                        beta[idx,i])
    write.table(select.cau,file = paste0(cur.dir,eth[i],"/select.cau_rho",l,"_",i1),row.names = F,col.names = F,quote=F)
  }
  
}








#all SNPs causal, genetic correlation 0.8, alpha = 0
i1 = 4
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
for(l in 1:3){
  eth <- c("EUR","AFR","AMR","EAS","SAS")
  cau.snp.infor <- cau.snp.infor.list[[l]]
  library(data.table)
  #generate the effect-size for all the SNPs
  #even if the SNPs have 0 MAF in a particular population the effectsize will aslo be generated
  #it's okay since only the extracted genotype are those only existing in a particular population
  
  #get the number of causal SNPs in EUR
  
  n.total.snp <- nrow(cau.snp.infor)
  
  
  sigma = 0.4
  gr12 = 0.8
  gr13 = 0.8
  gr14 = 0.8
  gr15 = 0.8
  gr23 = 0.8
  gr24 = 0.8
  gr25 = 0.8
  gr34 = 0.8
  gr35 = 0.8
  gr45 = 0.8
  n = as.numeric(n.total.snp)
  n1 = n
  n2 = n
  n3 = n
  n4 = n
  n5 = n
  Sigma <- GenSigma(sigma,n1,n2,n3,n4,n5,
                    gr12,gr13,gr14,gr15,
                    gr23,gr24,gr25,
                    gr34,gr35,
                    gr45)
  save(Sigma,file = paste0(cur.dir,"/causal_Sigma_",i1,".rdata"))
  library(mvtnorm)
  #gamma represent generate scale
  #u represent standarized scale effect-size
  #GCTA needs u
  set.seed(666)
  gamma =  rmvnorm(n.total.snp,c(0,0,0,0,0),
                  sigma=Sigma)
  colnames(gamma) <- paste0("gamma_",c("EUR","AFR","AMR","EAS","SAS"))
  #since this is the effect-size under standardized scale
  #to use GCTA, we need the original scale effect size
  EUR.bi <- cau.snp.infor$EUR>=0.01&
    cau.snp.infor$EUR<=0.99
  AFR.bi <-  cau.snp.infor$AFR>=0.01&
    cau.snp.infor$AFR<=0.99
  AMR.bi <-  cau.snp.infor$AMR>=0.01&
    cau.snp.infor$AMR<=0.99
  EAS.bi <-  cau.snp.infor$EAS>=0.01&
    cau.snp.infor$EAS<=0.99
  SAS.bi <-  cau.snp.infor$SAS>=0.01&
    cau.snp.infor$SAS<=0.99
  eth.bi <- cbind(EUR.bi,AFR.bi,AMR.bi,EAS.bi,SAS.bi)
  library(dplyr)
  MAF <- cau.snp.infor %>% select(EUR,AFR,AMR,EAS,SAS)
  MAF <- as.data.frame(MAF)
  for(i in 1:5){
    idx <- which(eth.bi[,i]==1)
    
    u = gamma[,i]*sqrt(2*MAF[,i]*(1-MAF[,i]))    
    #scale factor for heritability
    total.heritabillity = sum(u^2)
    
    u = u*sqrt(0.4/total.heritabillity)
    
    #beta.ori <- beta[idx,i]/sqrt(2*MAF[idx,i]*(1-MAF[idx,i]))
    #scale beta back to standarized scale
    # u =  beta_scale[idx]*sqrt(2*MAF[idx,i]*(1-MAF[idx,i]))
        
    
    select.cau <- cbind(cau.snp.infor[idx,1],
                        u[idx])
    write.table(select.cau,file = paste0(cur.dir,eth[i],"/select.cau_rho",l,"_",i1),row.names = F,col.names = F,quote=F)
  }
  
}


#all SNPs causal, genetic correlation 0.8, alpha = -0.25
i1 = 5
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
for(l in 1:3){
  eth <- c("EUR","AFR","AMR","EAS","SAS")
  cau.snp.infor <- cau.snp.infor.list[[l]]
  library(data.table)
  #generate the effect-size for all the SNPs
  #even if the SNPs have 0 MAF in a particular population the effectsize will aslo be generated
  #it's okay since only the extracted genotype are those only existing in a particular population
  
  #get the number of causal SNPs in EUR
  
  n.total.snp <- nrow(cau.snp.infor)
  
  
  sigma = 0.4
  gr12 = 0.8
  gr13 = 0.8
  gr14 = 0.8
  gr15 = 0.8
  gr23 = 0.8
  gr24 = 0.8
  gr25 = 0.8
  gr34 = 0.8
  gr35 = 0.8
  gr45 = 0.8
  n = as.numeric(n.total.snp)
  n1 = n
  n2 = n
  n3 = n
  n4 = n
  n5 = n
  Sigma <- GenSigma(sigma,n1,n2,n3,n4,n5,
                    gr12,gr13,gr14,gr15,
                    gr23,gr24,gr25,
                    gr34,gr35,
                    gr45)
  save(Sigma,file = paste0(cur.dir,"/causal_Sigma_",i1,".rdata"))
  library(mvtnorm)
  #beta represent effect under alpha = -0.25 scale
  #u represent standarized scale effect-size
  #GCTA needs u
  set.seed(666)
  gamma =  rmvnorm(n.total.snp,c(0,0,0,0,0),
                  sigma=Sigma)
  colnames(gamma) <- paste0("gamma_",c("EUR","AFR","AMR","EAS","SAS"))
  #since this is the effect-size under standardized scale
  #to use GCTA, we need the original scale effect size
  EUR.bi <- cau.snp.infor$EUR>=0.01&
    cau.snp.infor$EUR<=0.99
  AFR.bi <-  cau.snp.infor$AFR>=0.01&
    cau.snp.infor$AFR<=0.99
  AMR.bi <-  cau.snp.infor$AMR>=0.01&
    cau.snp.infor$AMR<=0.99
  EAS.bi <-  cau.snp.infor$EAS>=0.01&
    cau.snp.infor$EAS<=0.99
  SAS.bi <-  cau.snp.infor$SAS>=0.01&
    cau.snp.infor$SAS<=0.99
  eth.bi <- cbind(EUR.bi,AFR.bi,AMR.bi,EAS.bi,SAS.bi)
  library(dplyr)
  MAF <- cau.snp.infor %>% select(EUR,AFR,AMR,EAS,SAS)
  MAF <- as.data.frame(MAF)
  for(i in 1:5){
    idx <- which(eth.bi[,i]==1)
    
    #scale factor for heritability
    u = gamma[,i]*(2*MAF[,i]*(1-MAF[,i]))^0.375
    #scale factor for heritability
    total.heritabillity = sum(u^2)
    
    u = u*sqrt(0.4/total.heritabillity)
    
    #beta.ori <- beta[idx,i]/sqrt(2*MAF[idx,i]*(1-MAF[idx,i]))
    #scale beta back to standarized scale
    # u =  beta_scale[idx]*sqrt(2*MAF[idx,i]*(1-MAF[idx,i]))
    
    
    select.cau <- cbind(cau.snp.infor[idx,1],
                        u[idx])
    
    write.table(select.cau,file = paste0(cur.dir,eth[i],"/select.cau_rho",l,"_",i1),row.names = F,col.names = F,quote=F)
  }
  
}










