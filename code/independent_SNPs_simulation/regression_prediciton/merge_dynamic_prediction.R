#goal: merge dynamic prediction
#original coeficients
i3 = 1
pop.ind = 2
i1 = 1
n.s <- 100
r2.vad.mat <- matrix(0,n.s,4)
n.snp.mat <- matrix(0,n.s,4)
prop.mat <- matrix(0,n.s,4)

p.thr.mat <- matrix(0,n.s,4)
p.thr2.mat <- matrix(0,n.s,4)
r2.test.mat <- matrix(0,n.s,4)
pthr.train <- result[,1]
pthr.ref <- result[,2]
p.thr <- c(1E-8,5E-8,1E-7,5E-7,1E-6,5E-6,1E-5,5E-5,1E-4,5E-4,1E-3,5E-3,1E-2,0.1,0.3,0.5)
p.thr2 <- p.thr
r2.all.test <- matrix(0,length(p.thr)*length(p.thr2),4*n.s)
#i3 is method i3=1 represent orignal coef
#i3 =2 represent normal prior coef
#i1 index of simulation
#pop.ind 2 is African
#pop.ind 3 is LAT
temp <- 1
temp2 = 1 
gr <- 2
for(i3 in 1:2){
  for(pop.ind in 2:3){
  
    for(i1 in 1:100){
      load(paste0("./multi_ethnic/result/Dy_result_",i1,"_",pop.ind,"_",i3,"_",gr))
      idx <- which.max(result[,3])
      p.thr.mat[i1,temp] <- result[idx,1]
      p.thr2.mat[i1,temp] <- result[idx,2]
      r2.test.mat[i1,temp] <- result[idx,3]
      r2.vad.mat[i1,temp] <- result[idx,4]
      n.snp.mat[i1,temp] <- result[idx,5]
      prop.mat[i1,temp] <- result[idx,6]
      r2.all.test[,temp2] <- result[,3]
      temp2= temp2+1
    } 
    temp <- temp+1 
     }
 
}
r2.vad.mat
r2.test.mat


r2.plot.test <- matrix(0,length(p.thr)*length(p.thr2),4)
for(i in 1:4){
  temp <- r2.all.test[,c(1:n.s)+n.s*(i-1)]
  r2.plot.test[,i] <- rowMeans(temp)
}

data <- cbind(result[,c(1,2)],r2.plot.test)
colnames(data) <- c("pthr_train",
                    "pthr_ref",
                    "r2_AFR_ori",
                    "r2_AFR_EB",
                    "r2_LAC_ori",
                    "r2_LAC_EB")
save(data,file =paste0("./multi_ethnic/result/test_result_ori_EB.Rdata"))

colMeans(r2.vad.mat)
colMeans(n.snp.mat)
colMeans(prop.mat)


#
Pfunction <- fucntion(beta, beta_low,beta_high){
  sd_beta = (log(beta_high)-log(beta_low))/(2*qnorm(0.975))
  z = log(beta)/sd_beta
  p = 2*pnorm(-abs(z))
}
#calibrated coefficents
i3 = 2



if(i3==1){
  result <- LDPDy(y_all,
                  beta.train,
                  p.train,
                  p.thr,
                  pop.ind,
                  beta.ref,
                  p.ref,
                  alpha)
  save(result,file=paste0("./multi_ethnic/result/Dy_result_",i1,"_",pop.ind,"_",i3))  
  
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
                   alpha)
  save(result,file=paste0("./multi_ethnic/result/Dy_result_",i1,"_",pop.ind,"_",i3))  
}            

