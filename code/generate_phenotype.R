#Goal: Generate phenotypes data for
#simulate phenotypes data for AFR, EUR, LAC
#MAF based on 1000KG
#sample size EUR n =120000
#sample size AFR n = 18000
#sample size LAC n = 18000
#heritability for EUR 0.5
#heritability for AFR 0.5
#heritability for LAC 0.5
#5000 causal SNPs for each population
#4000 shared causal SNPs
#Genetic correlation for the between EUR and AFR is 0.4
#Genetic correlation for the between EUR and LAC is 0.6
#Genetic correlation for the between LAC and AFR is 0.6
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
#i2 is for genetic correlation structure
i2 = as.numeric(args[[2]])
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis')
load(paste0("./multi_ethnic/result/pruned_geno/geno_",1))

GenSigma <- function(her1,her2,her3,
                     gr12,gr13,gr23){
  vecsigma <- c(her1,
                gr12*sqrt(her1*her2),
                gr13*sqrt(her1*her3),
                gr12*sqrt(her1*her2),
                her2,
                gr23*sqrt(her2*her3),
                gr13*sqrt(her1*her3),
                gr23*sqrt(her2*her3),
                her3)
  Sigma <- matrix(vecsigma,3,3)
  return(Sigma)
          
}

gr.set <- list(c(0.4,0.6,0.6),
               c(0.6,0.4,0.6))

n.shared <- 4000
n.nonshare <- 1000
n.cau <- n.shared+n.nonshare
her.EUR <- 0.5
her.AFR <- 0.5
her.LAC <- 0.5
gr.vec <- gr.set[[i2]]
gr.EUR.AFR <- gr.vec[1]
gr.EUR.LAC <- gr.vec[2]
gr.LAC.AFR <- gr.vec[3]
Sigma <- GenSigma(her.EUR/n.cau,
                  her.AFR/n.cau,
                  her.LAC/n.cau,
                  gr.EUR.AFR,
                  gr.EUR.LAC,
                  gr.LAC.AFR)
library(mvtnorm)

#since the SNP MAF data is permutated
#for EUR and AFR, we will use 1:4000 as shared SNPs
#use 4000:5000 as nonshared SNPs for EUR
#5000:6000 as nonshared SNPs for AFR
#6000:7000 as nonshared SNPs for LAC


#shared beta coefcients
#######################################
# set.seed(666)
# beta_shared <- rmvnorm(n.shared,c(0,0,0),
#                         sigma=Sigma)
# # #non shared beta coefcients
#  all_nonshare <- 3*n.nonshare
#  beta_nonshared <- matrix(0,all_nonshare,3)
#  beta_nonshared[1:n.nonshare,1] <- rnorm(n.nonshare,0,sd=sqrt(her.EUR/n.cau))
#  beta_nonshared[(n.nonshare+1):(n.nonshare*2),2] <- rnorm(n.nonshare,0,sd=sqrt(her.AFR/n.cau))
#  beta_nonshared[(n.nonshare*2+1):(all_nonshare),3] <- rnorm(n.nonshare,0,sd=sqrt(her.LAC/n.cau))
#  beta <- rbind(beta_shared,beta_nonshared)
# # #save effect size for the causal SNPs
#  save(beta,file = paste0("./multi_ethnic/result/pruned_geno/effect_size_",i2,".Rdata"))
############################################



#beta is standardized genotype effect size
#b is orignal effect size
load("/spin1/users/zhangh24/KG.vcf/MAF_result/pruned_MAF_permu.Rdata")
MAF.EUR <- pruned.snp.permu$MAF.EUR
MAF.AFR <- pruned.snp.permu$MAF.AFR
MAF.LAC <- pruned.snp.permu$MAF.LAC
all_nonshare <- 3*n.nonshare
#n.all are all the SNPs in the three populations with non-zero effects
n.all <- n.shared+all_nonshare
load(paste0("./multi_ethnic/result/pruned_geno/effect_size_",i2,".Rdata"))
b <- beta
b[,1] <- beta[,1]/sqrt(2*(MAF.EUR*(1-MAF.EUR))[1:n.all])
b[,2] <- beta[,2]/sqrt(2*(MAF.AFR*(1-MAF.AFR))[1:n.all])
b[,3] <- beta[,3]/sqrt(2*(MAF.LAC*(1-MAF.LAC))[1:n.all])






set.seed(i1)
    n.EUR <- nrow(genotype[[1]])
    n.AFR <- nrow(genotype[[2]])
    n.LAC <- nrow(genotype[[3]])
    
    y_EUR <- genotype[[1]]%*%b[,1]+rnorm(n.EUR,0,sd=sqrt(1-her.EUR))
    y_AFR <- genotype[[2]]%*%b[,2]+rnorm(n.AFR,0,sd=sqrt(1-her.AFR))
    y_LAC <- genotype[[3]]%*%b[,3]+rnorm(n.LAC,0,sd=sqrt(1-her.LAC))
    
    y <- list(y_EUR,y_AFR,y_LAC)
    save(y,file = paste0("./multi_ethnic/result/y_",i1,"_",i2))
    
 
    
    


# try <-  genotype_EUR%*%b[1:5000,1]
# mean(try)




