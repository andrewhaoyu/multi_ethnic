#goal: implement the dynamic p-value methods
#Postbeta do EB on the standarized scale
#then transform back to orignal scale uses MAF
PostBeta <- function(beta,Sigma,Sigma0,MAF.train,
                     MAF.ref){
  n <- length(beta)
  beta_post <- solve(solve(Sigma)+solve(Sigma0))%*%(solve(Sigma)%*%beta)
  beta_post[1] <- beta_post[1]/sqrt(2*MAF.train*(1-MAF.train))
  beta_post[2] <- beta_post[2]/sqrt(2*MAF.ref*(1-MAF.ref))
  return(beta_post)
}

#p.thr for train
#p.thr2 for ref
LDPDyW <- function(y_all,
                  beta.train,
                  sd.train,
                  p.train,
                  p.thr,
                  pop.ind,
                  beta.ref,
                  sd.ref,
                  p.ref,
                  p.thr2,
                  MAF.train,
                  MAF.ref
                  ){
  y = y_all[[pop.ind]]
  n.sub <- length(y)
  n.train <- n.sub*10/12
  n.test <- n.sub/12
  n.vad <- n.sub/12
  y.test <- y[n.train+(1:n.test)]
  y.vad <- y[(n.train+n.test)+(1:n.vad)]
  n.cat <- length(p.thr2)*length(p.thr)
  r2.vad <- r2.test <- rep(0,n.cat)
  
  prop <- rep(0,n.cat)
  n.snp.sec <- rep(0,n.cat)
  pthr2.vec <- rep(0,n.cat)
  pthr.vec <- rep(0,n.cat)
  prs.mat <- matrix(0,n.test+n.vad,n.cat)

  #####create posterior beta for different calibration point
  beta.train.update<- matrix(0,length(beta.train),
                              n.cat)
  for(l in 1:n.cat){
    beta.train.update[,l] <- beta.train
  }
  beta.all <- cbind(beta.train,beta.ref)
  #beta standarized score
  beta.train.sd <- beta.train*sqrt(2*MAF.train*(1-MAF.train))
  beta.ref.sd <- beta.ref*sqrt(2*MAF.ref*(1-MAF.ref))
  
  sd.all <- cbind(sd.train,sd.ref)
  sd.train.sd <- sd.train*sqrt(2*MAF.train*(1-MAF.train))
  sd.ref.sd <- sd.ref*sqrt(2*MAF.ref*(1-MAF.ref))
  
  MAF.all <- cbind(MAF.train,MAF.ref)
  beta.all.sd <- cbind(beta.train.sd,
                       beta.ref.sd)
  sd.all.sd <- cbind(sd.train.sd,
                     sd.ref.sd)
  temp <- 1
  for(k in 1:length(p.thr2)){
    print(k)
    for(j in 1:length(p.thr)){
      print(j)
     #find the corresponding SNPs
      idx <- which(p.train<=p.thr[j]|
                     p.ref<=p.thr2[k])
      #if less than 3 snp, can't calculate covariance matrix
      if(length(idx)>=3)
      {
        Sigma0 <- cov(beta.all.sd[idx,])
        for(i in 1:length(idx)){
          Sigma <- diag(sd.all.sd[idx[i],]^2)
          beta.train.update[idx[i],temp]  <- PostBeta(beta.all.sd[idx[i],],
                                                      Sigma,
                                                      Sigma0,
                                                      MAF.train[idx[i]],
                                                      MAF.ref[idx[i]]
          )[1]
        }
        
      }
      
      temp <- temp+1
    }
    
  }
  
  
  n.snp <- 0
  
  #p.update <- apply(cbind(p.train*alpha[k],p.ref),1,min)
  for(i in 1:500){
    print(i)
    load(paste0("./multi_ethnic/result/pruned_geno/geno_",i))
    file.snp <- ncol(genotype[[pop.ind]])
    for(k in 1:length(p.thr2)){
      for(j in 1:file.snp){
        #the p-value should be smaller than either of the two
        if((p.ref[n.snp+j]>=p.thr2[k])){
          jdx <- which((p.train[n.snp+j]<=p.thr))  
        }else{
          jdx <- c(1:length(p.thr))
        }
        
        if(length(jdx)!=0){
          #different p.value thr (jdx) get different beta
          #different alpha (k) get differnt beta
          prs.temp <- matrix(0,n.test+n.vad,length(jdx))
          for(l in 1:length(jdx)){
            prs.temp[,l] <- genotype[[pop.ind]][n.train+(1:(n.test+n.vad)),j]*beta.train.update[n.snp+j,jdx+(k-1)*length(p.thr),drop=F][,l]  
          }
          
            prs.mat[,jdx+(k-1)*length(p.thr)] <- prs.mat[,jdx+(k-1)*length(p.thr)]+
            prs.temp
          }
      }
      }
    n.snp <- n.snp+file.snp
  }
  
  ###calculate the proportion of causal SNPs
  ###j is the reference population
  ###k is the target population
  temp <- 1
  for(j in 1:length(p.thr2)){
    for(k in 1:length(p.thr)){
      idx <- which((p.train<=p.thr[k])|
                     (p.ref<=p.thr2[j]))
      if(length(idx)==0){
        n.snp.sec[temp] <- 0
        prop[temp] <- 0
        r2.test[temp] <- 0
        r2.vad[temp] <- 0
        pthr2.vec[temp] <- p.thr2[j]
        pthr.vec[temp] <- p.thr[k]
        temp <- temp+1
      }else{
        n.snp.sec[temp] <- length(idx)
        cau.snp <- c(1:4000,4000+c(1:1000)+1000*(pop.ind-1))
        prop[temp] <- sum(idx%in%cau.snp)/length(idx)
        prs.test <- prs.mat[(1:n.test),temp]
        prs.vad <- prs.mat[n.test+(1:n.vad),temp]
        model1 <- lm(y.test~prs.test)
        r2.test[temp] <- summary(model1)$adj.r.squared
        model2 <- lm(y.vad~prs.vad)
        r2.vad[temp] <- summary(model2)$adj.r.squared
        pthr2.vec[temp] <- p.thr2[j]
        pthr.vec[temp] <- p.thr[k]
        temp <- temp+1
      }
    }
  }

  result <- data.frame(pthr.vec,
                       pthr2.vec,
                       r2.test,
                       r2.vad,
                       n.snp.sec,
                       prop)
  
  return(result)
}  






LDPDy <- function(y_all,
                  beta.train,
                  p.train,
                  p.thr,
                  pop.ind,
                  beta.ref,
                  p.ref,
                  p.thr2){
  y = y_all[[pop.ind]]
  n.sub <- length(y)
  n.train <- n.sub*10/12
  n.test <- n.sub/12
  n.vad <- n.sub/12
  y.test <- y[n.train+(1:n.test)]
  y.vad <- y[(n.train+n.test)+(1:n.vad)]
  n.cat <- length(p.thr2)*length(p.thr)
  r2.vad <- r2.test <- rep(0,n.cat)
  
  prop <- rep(0,n.cat)
  n.snp.sec <- rep(0,n.cat)
  pthr.vec <- rep(0,n.cat)
  pthr2.vec <- rep(0,n.cat)
  prs.mat <- matrix(0,n.test+n.vad,n.cat)
  
  ######create a p value matrix to contain all different situation
  n.snp <- 0
  
  
  #p.update <- apply(cbind(p.train*alpha[k],p.ref),1,min)
  for(i in 1:500){
    print(i)
    load(paste0("./multi_ethnic/result/pruned_geno/geno_",i))
    file.snp <- ncol(genotype[[pop.ind]])
    for(k in 1:length(p.thr2)){
      #p.update <- p.mat[,k]
      for(j in 1:file.snp){
        if((p.ref[n.snp+j]>=p.thr2[k])){
          jdx <- which((p.train[n.snp+j]<=p.thr))  
        }else{
          jdx <- c(1:length(p.thr))
        }
        
        #jdx <- which(p.update[n.snp+j]<=p.thr)
        if(length(jdx)!=0){
          prs.temp <- genotype[[pop.ind]][,j]*beta.train[n.snp+j]
          prs.mat[,jdx+(k-1)*length(p.thr)] <- prs.mat[,jdx+(k-1)*length(p.thr)]+
            prs.temp[n.train+(1:(n.test+n.vad))]
        }
      }
    }
    n.snp <- n.snp+file.snp
  }
  ###calculate the proportion of causal SNPs
  ###j is the reference population
  ###k is the target population
  temp <- 1
  for(j in 1:length(p.thr2)){
    for(k in 1:length(p.thr)){
      idx <- which((p.train<=p.thr[k])|
                     (p.ref<=p.thr2[j]))
      if(length(idx)==0){
        n.snp.sec[temp] <- 0
        prop[temp] <- 0
        r2.test[temp] <- 0
        r2.vad[temp] <- 0
       pthr2.vec[temp] <- p.thr2[j]
        pthr.vec[temp] <- p.thr[k]
        temp <- temp+1
      }else{
        n.snp.sec[temp] <- length(idx)
        cau.snp <- c(1:4000,4000+c(1:1000)+1000*(pop.ind-1))
        prop[temp] <- sum(idx%in%cau.snp)/length(idx)
        prs.test <- prs.mat[(1:n.test),temp]
        prs.vad <- prs.mat[n.test+(1:n.vad),temp]
        model1 <- lm(y.test~prs.test)
        r2.test[temp] <- summary(model1)$adj.r.squared
        model2 <- lm(y.vad~prs.vad)
        r2.vad[temp] <- summary(model2)$adj.r.squared
        pthr2.vec[temp] <- p.thr2[j]
        pthr.vec[temp] <- p.thr[k]
        temp <- temp+1
      }
    }
  }
  
  result <- data.frame(pthr.vec,
                       pthr2.vec,
                       r2.test,
                       r2.vad,
                       n.snp.sec,
                       prop)
  
  return(result)
}  


arg <- commandArgs(trailingOnly=T)
#i1 represent phenotype file
i1 <- as.numeric(arg[[1]])
#pop.ind represent population
pop.ind <- as.numeric(arg[[2]])
#i3 represent orginal method or EB
i3 <- as.numeric(arg[[3]])
#i4 represent genetic correlation pattern
gr <- as.numeric(arg[[4]])

setwd('/spin1/users/zhangh24/breast_cancer_data_analysis')
load("/spin1/users/zhangh24/KG.vcf/MAF_result/pruned_MAF_permu.Rdata")
MAF.EUR <- pruned.snp.permu$MAF.EUR
MAF.AFR <- pruned.snp.permu$MAF.AFR
MAF.LAC <- pruned.snp.permu$MAF.LAC
load(paste0("./multi_ethnic/result/pruned_geno/beta_all_",i1,"_",gr,".Rdata"))
load(paste0("./multi_ethnic/result/y_",i1,"_",gr))
y_all <- y
beta_train <- beta_result[[1]]
sd_train1 <- beta_train[,2]
p_train1 <- beta_train[,3]
sd_train2 <- beta_train[,5]
p_train2 <- beta_train[,6]
sd_train3 <- beta_train[8]
p_train3 <- beta_train[,9]
beta_train1 <- beta_train[,1]
beta_train2 <- beta_train[,4]
beta_train3 <- beta_train[,7]

p.thr <- c(10^-8,5E-8,10^-7,5E-7,10^-6,5E-6,10^-5,5E-5,10^-4,5E-4,10^-3,5E-3,10^-2,0.1,0.3,0.5)

p.thr2 <- c(10^-8,5E-8,10^-7,5E-7,10^-6,5E-6,10^-5,5E-5,10^-4,5E-4,10^-3,5E-3,10^-2,0.1,0.3,0.5)

beta.ref <- beta_train1
p.ref <- p_train1
sd.ref <- sd_train1
p.train <- beta_train[,3*pop.ind]
beta.train <- beta_train[,3*pop.ind-2]
sd.train <- beta_train[,3*pop.ind-1]
MAF.ref <- MAF.EUR

if(pop.ind==2){
  MAF.train = MAF.AFR
}else{
  MAF.train = MAF.LAC
}

if(i3==1){
  result <- LDPDy(y_all,
        beta.train,
        p.train,
        p.thr,
        pop.ind,
        beta.ref,
        p.ref,
        p.thr2)
save(result,file=paste0("./multi_ethnic/result/Dy_result_",i1,"_",pop.ind,"_",i3,"_",gr))  
  
}else{
  result <- LDPDyW(y_all,
                     beta.train,
                     sd.train,
                     p.train,
                     p.thr,
                     pop.ind,
                     beta.ref,
                     sd.ref,
                     p.ref,
                      p.thr2,
                      MAF.train,
                      MAF.ref)
  save(result,file=paste0("./multi_ethnic/result/Dy_result_",i1,"_",pop.ind,"_",i3,"_",gr))  
}            

