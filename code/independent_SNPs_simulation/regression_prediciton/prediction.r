#Goal: implement three prediction methods
#method1: LD pruning and threshold
#method2: weighted combination of the three
#method3: E-Bayes
arg <- commandArgs(trailingOnly=T)
#i1 represent phenotype file
i1 <- as.numeric(arg[[1]])
#pop.ind represent population type
pop.ind <- as.numeric(arg[[2]])
#gr represent gentic correlation
gr <- as.numeric(arg[[3]])
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis')
load(paste0("./multi_ethnic/result/y_",i1,"_",gr))
load(paste0("./multi_ethnic/result/pruned_geno/beta_all_",i1,"_",gr,".Rdata"))


y_all = y


#find all genotypes data files cutpoint
# all.cut <- rep(0,500)
# n.snp <- 0
# for(i1 in 1:500){
#   if(i1%%50==0){
#     print(i1)  
#   }
#   load(paste0("./multi_ethnic/result/pruned_geno/beta_estimate_",i1,".Rdata"))
#  temp <- nrow(beta_result[[1]])
#   n.snp <- n.snp +temp
#   all.cut[i1] <- n.snp
# }




#pop.ind population indicator 1 EUR, 2 AFR, 3 LAC
LDP <- function(y_all,beta.train,p.train,
                p.thr,pop.ind){
  y = y_all[[pop.ind]]
  n.sub <- length(y)
  n.train <- n.sub*10/12
  n.test <- n.sub/12
  n.vad <- n.sub/12
  y.test <- y[n.train+(1:n.test)]
  y.vad <- y[(n.train+n.test)+(1:n.vad)]
  r2.vad <- r2.test <- rep(0,length(p.thr))
  prop <- rep(0,length(p.thr))
  n.snp.sec <- rep(0,length(p.thr))
  prs.mat <- matrix(0,n.test+n.vad,length(p.thr))
  
  
  n.snp <- 0
  for(i in 1:500){
    print(i)
    load(paste0("./multi_ethnic/result/pruned_geno/geno_",i))
    file.snp <- ncol(genotype[[pop.ind]])
    for(j in 1:file.snp){
      jdx <- which(p.train[n.snp+j]<=p.thr)
        if(length(jdx)!=0){
          prs.temp <- genotype[[pop.ind]][,j]*beta.train[n.snp+j]
          prs.mat[,jdx] <- prs.mat[,jdx]+
          prs.temp[n.train+(1:(n.test+n.vad))]
        }
    }
    n.snp <- n.snp+file.snp
  }

  for(k in 1:length(p.thr)){
    idx <- which(p.train<=p.thr[k])
    if(length(idx)==0){
      n.snp.sec[k] <- 0
      prop[k] <- 0
      r2.test[k] <- 0
      r2.vad[k] <- 0
    }else{
      n.snp.sec[k] <- length(idx)
      cau.snp <- c(1:4000,4000+c(1:1000)+1000*(pop.ind-1))
      prop[k] <- sum(idx%in%cau.snp)/length(idx)
      prs.test <- prs.mat[(1:n.test),k]
      prs.vad <- prs.mat[n.test+(1:n.vad),k]
      model1 <- lm(y.test~prs.test)
      r2.test[k] <- summary(model1)$adj.r.squared
      model2 <- lm(y.vad~prs.vad)
      r2.vad[k] <- summary(model2)$adj.r.squared
    }
  }

  result <- list(n.snp.sec,prop,
                 r2.test,r2.vad,
                 prs.mat)
  return(result)
}  
#LDP EUR
LDPEUR <- function(y_all,
                beta.train,
                p.train,
                p.thr,
                pop.ind){
  y = y_all[[pop.ind]]
  n.sub <- length(y)
  n.train <- n.sub*10/12
  n.test <- n.sub/12
  n.vad <- n.sub/12
  y.test <- y[n.train+(1:n.test)]
  y.vad <- y[(n.train+n.test)+(1:n.vad)]
  r2.vad <- r2.test <- rep(0,length(p.thr))
  prop <- rep(0,length(p.thr))
  n.snp.sec <- rep(0,length(p.thr))
  prs.mat <- matrix(0,n.test+n.vad,length(p.thr))
  n.sub2 <- length(y_all[[2]])
  n.train2 <- n.sub2*10/12
  n.test2 <- n.sub2/12
  n.vad2 <- n.sub2/12
  prs.mat2 <- matrix(0,n.test2+n.vad2,length(p.thr))
  n.sub3 <- length(y_all[[3]])
  n.train3 <- n.sub3*10/12
  n.test3 <- n.sub3/12
  n.vad3 <- n.sub3/12
  prs.mat3 <- matrix(0,n.test2+n.vad2,length(p.thr))
  
  n.snp <- 0
  for(i in 1:500){
    print(i)
    load(paste0("./multi_ethnic/result/pruned_geno/geno_",i))
    file.snp <- ncol(genotype[[pop.ind]])
    for(j in 1:file.snp){
      jdx <- which(p.train[n.snp+j]<=p.thr)
      if(length(jdx)!=0){
        prs.temp <- genotype[[pop.ind]][,j]*beta.train[n.snp+j]
        prs.mat[,jdx] <- prs.mat[,jdx]+
          prs.temp[n.train+(1:(n.test+n.vad))]
        prs.temp2 <- genotype[[2]][,j]*beta.train[n.snp+j] 
        prs.mat2[,jdx] <- prs.mat2[,jdx]+
          prs.temp2[n.train2+(1:(n.test2+n.vad2))]
        prs.temp3 <- genotype[[3]][,j]*beta.train[n.snp+j] 
        prs.mat3[,jdx] <- prs.mat3[,jdx]+
          prs.temp3[n.train3+(1:(n.test3+n.vad3))]
      }
    }
    n.snp <- n.snp+file.snp
  }
  
  for(k in 1:length(p.thr)){
    idx <- which(p.train<=p.thr[k])
    if(length(idx)==0){
      n.snp.sec[k] <- 0
      prop[k] <- 0
      r2.test[k] <- 0
      r2.vad[k] <- 0
    }else{
      n.snp.sec[k] <- length(idx)
      cau.snp <- c(1:4000,4000+c(1:1000)+1000*(pop.ind-1))
      prop[k] <- sum(idx%in%cau.snp)/length(idx)
      prs.test <- prs.mat[(1:n.test),k]
      prs.vad <- prs.mat[n.test+(1:n.vad),k]
      model1 <- lm(y.test~prs.test)
      r2.test[k] <- summary(model1)$adj.r.squared
      model2 <- lm(y.vad~prs.vad)
      r2.vad[k] <- summary(model2)$adj.r.squared
    }
  }
  
  result <- list(n.snp.sec,prop,
                 r2.test,r2.vad,
                 prs.mat,
                 prs.mat2,
                 prs.mat3)
  return(result)
}  
p.thr <- c(10^-8,5E-8,10^-7,5E-7,10^-6,5E-6,10^-5,5E-5,10^-4,5E-4,10^-3,5E-3,10^-2,0.1,0.3,0.5)
beta.train <- beta_result[[1]][,3*pop.ind-2]
beta.test <- beta_result[[2]][,3*pop.ind-2]
beta.vad <- beta_result[[3]][,3*pop.ind-2]

p.train <- beta_result[[1]][,3*pop.ind]

if(pop.ind==1){
  LDP.result <-  LDPEUR(y_all,beta.train,p.train,p.thr,pop.ind)
}else{
  LDP.result <-  LDP(y_all,beta.train,p.train,p.thr,pop.ind)
}


save(LDP.result,file = paste0("./multi_ethnic/result/LDP.result_",i1,"_",pop.ind,"_",gr))
  

# 
#   
# LDP <- function(y_all,beta.train,beta.test,beta.vad,p.train,
#                 p.thr,pop.ind){
#   y = y_all[[pop.ind]]
#   n.sub <- length(y)
#   n.train <- n.sub*10/12
#   n.test <- n.sub/12
#   n.vad <- n.sub/12
#   y.test <- y[n.train+(1:n.test)]
#   y.vad <- y[(n.train+n.test)+(1:n.vad)]
#   r2.vad <- r2.test <- rep(0,length(p.thr))
#   prop <- rep(0,length(p.thr))
#   n.snp.sec <- rep(0,length(p.thr))
#   prs.mat <- matrix(0,n.test+n.vad,length(p.thr))
#   
#   for(k in 1:length(p.thr)){
#     #print(k)
#     idx <- which(p.train<=p.thr[k])
#     if(length(idx)==0){
#       n.snp.sec[k] <- 0
#       prop[k] <- 0
#       r2.test[k] <- 0
#       r2.vad[k] <- 0
#     }else{
#       n.snp.sec[k] <- length(idx)
#       cau.snp <- c(1:4000,4000+c(1:1000)+1000*(pop.ind-1))
#       prop[k] <- sum(idx%in%cau.snp)/length(idx)
#       prs <- rep(0,n.sub)
#       temp.result <- Findfilenum(idx,all.cut)
#       filenum <- temp.result[[1]]
#       colnum <- temp.result[[2]]
#       if(is.null(filenum)){
#         r2.test[k] <- NULL
#         r2.vad[k] <- NULL
#       }else{
#         
#         for(i in 1:length(filenum)){
#           print(i)
#           if(i==1){
#             #####load the first genotype file
#             file.temp <- filenum[i]
#             load(paste0("./multi_ethnic/result/pruned_geno/geno_",file.temp))
#           }else if(file.temp!=filenum[i]){
#             #avoid reload genotype data
#             load(paste0("./multi_ethnic/result/pruned_geno/geno_",filenum[i]))
#             file.temp <- filenum[i]
#           }
#           #take the corresponding column number for selected SNP
#           geno <- genotype[[pop.ind]][,colnum[i]]
#           prs.temp <- beta.train[idx[i]]*geno
#           prs <- prs+prs.temp  
#         }
#       }
#       prs.test <- prs[n.train+(1:n.test)]
#       prs.vad <- prs[n.train+n.test+(1:n.vad)]
#       model1 <- lm(y.test~prs.test)
#       r2.test[k] <- summary(model1)$adj.r.squared
#       model2 <- lm(y.vad~prs.vad)
#       r2.vad[k] <- summary(model2)$adj.r.squared
#       prs.mat[,k] <- c(prs.test,prs.vad)
#     }
#     
#   }
#   result <- list(n.snp.sec,prop,
#                  r2.test,r2.vad,
#                  prs.mat)
#   return(result)
#   
# }
# 

  


# #find the file number and column number for significant SNP
# Findfilenum <- function(idx, all.cut){
#   n.idx <- length(idx)
#   if(n.idx==0){
#     filenum = NULL
#     colnum = NULL
#   }else{
#     filenum <- rep(0,n.idx)
#     colnum <- rep(0,n.idx)
#     for(j in 1:length(idx)){
#       for(i in 1:500){
#         
#         if(idx[j] <=all.cut[i]){
#           filenum[j] <- i
#           if(i==1){
#             colnum[j] <- idx[j]
#           }else{
#             colnum[j] <- idx[j]-all.cut[i-1]
#           }
#           break;}
#         }
#     }
#   }
#   return(list(filenum,colnum))
# }


