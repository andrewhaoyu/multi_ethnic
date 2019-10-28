#Goal: merge prediction LDP result
#simulation number 100
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis')
n.s <- 100
p.thr <- c(10^-8,5E-8,10^-7,5E-7,10^-6,5E-6,10^-5,5E-5,10^-4,5E-4,10^-3,5E-3,10^-2,0.1,0.3,0.5)
n.p <- length(p.thr)
n.snp.mat <- matrix(0,n.s,n.p)
prop.mat <- matrix(0,n.s,n.p)
r2.test.mat <- matrix(0,n.s,n.p)
r2.vad.mat <- matrix(0,n.s,n.p)
#the best testing threshold adjusted r2 on validation data
best.vad <- rep(0,n.s)
gr <- 2
for(pop.ind in 1:3){
for(i1 in 1:n.s){
 load( paste0("./multi_ethnic/result/LDP.result_",i1,"_",pop.ind,"_",gr))
      #load(paste0("./multi_ethnic/result/LDP.result_",i1,"_",pop.ind))
      n.snp.mat[i1,] <- LDP.result[[1]]
      prop.mat[i1,] <- LDP.result[[2]]
      r2.test.mat[i1,] <- LDP.result[[3]]
      r2.vad.mat[i1,] <- LDP.result[[4]]
      idx <- which.max(LDP.result[[3]])
      best.vad[i1] <- LDP.result[[4]][idx]
      
      
}
  ##put all the rows with 0 SNPs selected as NA
  for(k in 1:ncol(prop.mat)){
    idx <- which(prop.mat[,k]==0)
    prop.mat[idx,k] <- NA
  }
  #prop is the proportion 
  n.snp <- colMeans(n.snp.mat)
  prop <- colMeans(prop.mat,na.rm=T)
  r2.test <- colMeans(r2.test.mat)
  r2.vad <- colMeans(r2.vad.mat)
    LDP.result <- list(n.snp,
                       prop,
                       r2.test,
                       r2.vad,
                       best.vad)
    save(LDP.result, file = paste0("./multi_ethnic/result/LDP_summary_",pop.ind,"_",gr))
}




