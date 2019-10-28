#goal: weighted combination of PRS methods
n.s <- 100
#vadlidation results for AFR and LAT
r2.2 <- rep(0,n.s)
r2.3 <- rep(0,n.s)
gr <- 2
for(i1 in 1:100){
  setwd('/spin1/users/zhangh24/breast_cancer_data_analysis')
  pop.ind <- 1
  load(paste0("./multi_ethnic/result/y_",i1,"_",gr))
  y_all = y
  load(paste0("./multi_ethnic/result/LDP.result_",i1,"_",pop.ind,"_",gr))
  prs.mat2 <- LDP.result[[6]]
  prs.mat3 <- LDP.result[[7]]
  r2.test <- LDP.result[[3]]
  idx.best <-which.max(r2.test)
  prs.EUR2 <- prs.mat2[,idx.best]
  prs.EUR3 <- prs.mat3[,idx.best]
  
  #weighted combination on AFR
  pop.ind <- 2
  load(paste0("./multi_ethnic/result/LDP.result_",i1,"_",pop.ind,"_",gr))
  r2.test <- LDP.result[[3]]
  idx.best <-which.max(r2.test)
  prs.AFR <- LDP.result[[5]][,idx.best]
  prs.AFR.test <- prs.AFR
  y = y_all[[pop.ind]]
  n.sub <- length(y)
  n.train <- n.sub*10/12
  n.test <- n.sub/12
  n.vad <- n.sub/12
  y.test <- y[n.train+(1:n.test)]
  prs.AFR.test <- prs.AFR[(1:n.test)]
  prs.EUR2.test <- prs.EUR2[(1:n.test)]
  model1 <- lm(y.test~prs.AFR.test+prs.EUR2.test)
  weight <- coef(model1)[2:3]
  y.vad <- y[n.train+n.test+(1:n.vad)]
  prs.AFR.vad <- prs.AFR[n.test+(1:n.vad)]
  prs.EUR2.vad <- prs.EUR2[n.test+(1:n.vad)]
  prs.weight <- cbind(prs.AFR.vad,prs.EUR2.vad)%*%weight
  model2 <- lm(y.vad~prs.weight)
  r2.raw <- summary(model2)$r.squared
  r2.adusted.2 <- 1-(1-r2.raw)*(n.vad-1)/(n.vad-2-1)
  r2.2[i1] <- r2.adusted.2
  #weighted combination on LAC
  pop.ind <- 3
  load(paste0("./multi_ethnic/result/LDP.result_",i1,"_",pop.ind,"_",gr)) 
  r2.test <- LDP.result[[3]]
  idx.best <-which.max(r2.test)
  prs.LAC <- LDP.result[[5]][,idx.best]
  prs.LAC.test <- prs.LAC
  y = y_all[[pop.ind]]
  n.sub <- length(y)
  n.train <- n.sub*10/12
  n.test <- n.sub/12
  n.vad <- n.sub/12
  y.test <- y[n.train+(1:n.test)]
  prs.LAC.test <- prs.LAC[(1:n.test)]
  prs.EUR3.test <- prs.EUR3[(1:n.test)]
  model1 <- lm(y.test~prs.LAC.test+prs.EUR3.test)
  weight <- coef(model1)[2:3]
  y.vad <- y[n.train+n.test+(1:n.vad)]
  prs.LAC.vad <- prs.LAC[n.test+(1:n.vad)]
  prs.EUR3.vad <- prs.EUR3[n.test+(1:n.vad)]
  prs.weight <- cbind(prs.LAC.vad,prs.EUR3.vad)%*%weight
  model2 <- lm(y.vad~prs.weight)
  r2.raw <- summary(model2)$r.squared
  r2.adusted.3 <- 1-(1-r2.raw)*(n.vad-1)/(n.vad-2-1)
  r2.3[i1] <- r2.adusted.3
  
  
}

