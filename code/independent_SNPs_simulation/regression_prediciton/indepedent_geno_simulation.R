#simulate phenotypes data for AFR and EUR
#MAF based on 1000KG
#sample size EUR n =100000
#sample size AFR n = 15000
#heritability for EUR 0.8
#heritability for AFR 0.8
#5000 causal SNPs for each population
#4000 shared causal SNPs
#Genetic correlation for the shared SNPs is 0.6
#1000 SNPs for each as independent causal
# load("/data/zhangh24/KG.vcf/MAF_result/pruned_MAF.Rdata")
#  set.seed(666)
#  n.snp <- nrow(pruned.snp.clean)
#  pruned.snp.permu <- pruned.snp.clean[sample(c(1:n.snp),n.snp),]
#  save(pruned.snp.permu,file = "/data/zhangh24/KG.vcf/MAF_result/pruned_MAF_permu.Rdata")

load("/data/zhangh24/KG.vcf/MAF_result/pruned_MAF_permu.Rdata")
set.seed(666)
idx <- which.min(pruned.snp.permu$all.snp.AFR.MAF.AFR[1:5000])
pruned.snp.permu[idx,]
library(RcppArmadillo)
MAF.EUR <- pruned.snp.permu$MAF.EUR
MAF.AFR <- pruned.snp.permu$all.snp.AFR.MAF.AFR
n.EUR <- 120000
n.AFR <- 18000
n.train.EUR <- 100000
n.test.EUR <- 10000
n.valid.EUR <- 10000
n.train.AFR <- 15000
n.test.AFR <- 1500
n.valid.AFR <- 1500
n.shared <- 4000
n.nonshare <- 1000
n.cau <- n.shared+n.nonshare
her.EUR <- 0.8
her.AFR <- 0.8
gr <- 0.6
library(mvtnorm)


Sigma <- matrix(c(her.EUR/n.cau,
                  gr*sqrt(her.EUR/n.cau*her.AFR/n.cau),
                  gr*sqrt(her.EUR/n.cau*her.AFR/n.cau),
                  her.AFR/n.cau),2,2)



#since the SNP MAF data is permutated
#for EUR and AFR, we will use 1:4000 as shared SNPs
#use 4000:5000 as nonshared SNPs for EUR
#5000:6000 as nonshared SNPs for AFR
G_EUR_shared <- matrix(rbinom(n.shared*n.EUR,2,MAF.EUR[1:n.shared]),n.EUR,n.shared,byrow = T)
G_AFR_shared <- matrix(rbinom(n.shared*n.AFR,2,MAF.AFR[1:n.shared]),n.AFR,n.shared,byrow = T)
#for easyness of coding, we will random simulate 2*nonshared SNPs
#but put second half for AFR as nonshare
two_nonshare = 2*n.nonshare
G_EUR_nonshare <- matrix(rbinom(two_nonshare*n.EUR,2,MAF.EUR[(n.shared+1):(n.shared+two_nonshare)]),n.EUR,two_nonshare,byrow = T)
G_EUR_effect <- cbind(G_EUR_shared,G_EUR_nonshare)
G_AFR_nonshare <- matrix(rbinom(two_nonshare*n.AFR,2,MAF.AFR[(n.shared+1):(n.shared+two_nonshare)]),n.AFR,two_nonshare,byrow = T)
G_AFR_effect <- cbind(G_AFR_shared,G_AFR_nonshare)

#shared beta coefcients
beta_shared <- rmvnorm(n.shared,c(0,0),
                sigma=Sigma)
#non shared beta coefcients
beta_nonshared <- matrix(0,two_nonshare,2)
beta_nonshared[1:n.nonshare,1] <- rnorm(n.nonshare,0,sd=sqrt(her.EUR/n.cau))
beta_nonshared[(n.nonshare+1):(two_nonshare),2] <- rnorm(n.nonshare,0,sd=sqrt(her.AFR/n.cau))
beta <- rbind(beta_shared,beta_nonshared)

y_EUR <- G_EUR_effect%*%beta[,1]+rnorm(n.EUR,0,sd=sqrt(1-her.EUR))
y_AFR <- G_AFR_effect%*%beta[,2]+rnorm(n.AFR,0,sd=sqrt(1-her.AFR))

y_EUR_train <- y_EUR[1:n.train.EUR]
y_AFR_train <- y_AFR[1:n.train.AFR]
y_EUR_test <- y_EUR[(n.train.EUR+1):(n.train.EUR+n.test.EUR)]
y_AFR_test <- y_AFR[(n.train.AFR+1):(n.train.AFR+n.test.AFR)]


####regression on the training and testing dataset
####record the summary level statistics
####first on the effect SNPs
####simulate the other SNPs and record the summary level statistics
n.snp <- nrow(pruned.snp.permu)
beta_summary_train <- matrix(0,n.snp,6)
beta_summary_test <- matrix(0,n.snp,6)
colnames(beta_summary_train) <- c("beta_EUR","sd_EUR","p_EUR",
                                  "beta_AFR","sd_AFR","p_AFR")
colnames(beta_summary_test) <- colnames(beta_summary_train)


FitLinearmodel <- function(y,x){
  model <- fastLm(X=cbind(1,x),y=y)
  if(is.na(coef(model)[2])){
    result <- c(0,1,1)
  }else{
    result <- coef(summary(model))[2,c(1,2,4)]  
  }
  
  return(result)
}
temp=1
for(i in 1:n.cau){
  if(i%%100==0){
    print(i)
  }
  beta_summary_train[temp,1:3] <- FitLinearmodel(y_EUR_train,G_EUR_effect[1:n.train.EUR,i])
  beta_summary_test[temp,1:3] <-  FitLinearmodel(y_EUR_test,G_EUR_effect[(n.train.EUR+1):(n.train.EUR+n.test.EUR),i])
  beta_summary_train[temp,4:6] <- FitLinearmodel(y_AFR_train,G_AFR_effect[1:n.train.AFR,i])
  beta_summary_test[temp,4:6] <-  FitLinearmodel(y_AFR_test,G_AFR_effect[(n.train.AFR+1):(n.train.AFR+n.test.AFR),i])
  temp = temp+1
}

for(i in 1:(n.snp-n.cau)){
  if(i%%100==0){
    print(i)
  }
  G_EUR <- rbinom(n.EUR,2,MAF.EUR[n.cau+i])
  G_AFR <- rbinom(n.AFR,2,MAF.AFR[n.cau+i])
  beta_summary_train[temp,1:3] <- FitLinearmodel(y_EUR_train,G_EUR[1:n.train.EUR])
  beta_summary_test[temp,1:3] <-  FitLinearmodel(y_EUR_test,G_EUR[(n.train.EUR+1):(n.train.EUR+n.test.EUR)])
  beta_summary_train[temp,4:6] <- FitLinearmodel(y_AFR_train,G_AFR[1:n.train.AFR])
  beta_summary_test[temp,4:6] <-  FitLinearmodel(y_AFR_test,G_AFR[(n.train.AFR+1):(n.train.AFR+n.test.AFR)])
  temp = temp+1
}


temp.result <- list(beta_summary_train,
                    beta_summary_test,
                     )


