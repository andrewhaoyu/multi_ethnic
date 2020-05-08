args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])

load("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/causal_genotype.rdata")
MAF <- rowSums(genotype)/(2*ncol(genotype))
genotype_standard <- genotype
for(i in 1:length(MAF)){
  genotype_standard[i,] <- (genotype[i,]-2*MAF[i])/(sqrt(2*MAF[i]*(1-MAF[i])))
}

genotype_standard_train = genotype_standard[,1:100000]
genotype_standard_test = genotype_standard[,100001:110000]

step = 100000/1000


result <- rep(0,step)
for(k in (1+(i1-1)*step):(i1*step)){
  print(k)
  result[k] <- max(cor(genotype_standard_train[,k],genotype_standard_test))  
}
save(result,file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/relatedness_",i1,".rdata"))


#merge the result
result_all <- rep(0,100000)
total <- 0
for(k in 1:1000){
  load(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/relatedness_",i1,".rdata"))
  temp <- length(result)
  result_all[total+(1:temp)] <- result
}

#temp <- apply(genotype_standard,1,sd)
#temp2 <- apply(genotype_standard,1,mean)

a <- matrix(rnorm(6),3,2)
cor(a)
