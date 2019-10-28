#Goal: merge summary level statistics
#get the total number of snps
n.snp <- 0
k <- 1
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis')
for(i1 in 1:500){
  if(i1%%50==0){
    print(i1)  
  }
  load(paste0("./multi_ethnic/result/pruned_geno/beta_estimate_",k,"_",i1,".Rdata"))
  temp <- nrow(beta_result[[1]])
  n.snp = n.snp +temp
}

for(k in 1:100){
  
  
  beta_summary_train <- matrix(0,n.snp,9)
  beta_summary_test <- matrix(0,n.snp,9)
  beta_summary_vad <- matrix(0,n.snp,9)
  
  
  colnames(beta_summary_train) <- c("beta_EUR","sd_EUR","p_EUR",
                                    "beta_AFR","sd_AFR","p_AFR", "beta_LAC","sd_LAC","p_LAC")
  n.snp <- 0
  #k represent phenotype
  #i1 represent genotype
  #gr represent genetic correlation
  gr = 2
  for(i1 in 1:500){
    if(i1%%50==0){
      print(i1)  
    }
    
    load(paste0("./multi_ethnic/result/pruned_geno/beta_estimate_",k,"_",i1,"_",gr,".Rdata"))
    temp <- nrow(beta_result[[1]])
    beta_summary_train[n.snp+(1:temp),] <- beta_result[[1]]
    beta_summary_test[n.snp+(1:temp),] <- beta_result[[2]]
    beta_summary_vad[n.snp+(1:temp),] <- beta_result[[3]]
    n.snp = n.snp +temp
  }
  
  beta_result <- list(beta_summary_train,
                      beta_summary_test,
                      beta_summary_vad)
  save(beta_result,file = paste0("./multi_ethnic/result/pruned_geno/beta_all_",k,"_",gr,".Rdata"))
  }



