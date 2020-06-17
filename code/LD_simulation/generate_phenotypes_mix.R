#generate phenotypes for different ethnic groups
#load the cuasal SNPs information
load("/data/zhangh24/multi_ethnic/result/LD_simulation/cau.snp.infor.rdata")
#order the causal SNPs by CHR and position

idx.order <- order(cau.snp.infor$CHR,
                   cau.snp.infor$position)
cau.snp.infor <- cau.snp.infor[idx.order,]
#head(cau.snp.infor)
#cau.snp.infor = cau.snp.infor[,-c(13,14,15,16)]
#head(cau.snp.infor)
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
n.total.snp <- nrow(cau.snp.infor)


GetEthnicProportion <- function(target,MAF){
  cau.snp.infor.shared=  cau.snp.infor%>% filter(
    get("EUR")>=MAF&
      get("EUR")<=(1-MAF)&
      get(target)>=MAF&
      get(target)<=1-MAF
  )
n.shared = nrow(cau.snp.infor.shared)  
cau.snp.infor.target=  cau.snp.infor%>% filter(
    get(target)>=MAF&
    get(target)<=1-MAF
)
n.target = nrow(cau.snp.infor.target)
return(1-n.shared/n.target)
}
eth <- c("AFR","AMR","EAS")
eth_proportion <- rep(0,length(eth))
for(k in 1:length(eth)){
  eth_proportion[k] = GetEthnicProportion(eth[k],0.05)
}


all.shared = Reduce(intersect,list(EUR.idx,EAS.idx,AMR.idx,AFR.idx))
n.all.shared = length(all.shared)

GenSigma <- function(sigma,n1,n2,n3,n4,
                     gr12,gr13,g14,
                     gr23,g24,g34){
  vecsigma <- sigma*c(1/n1,
                      gr12/sqrt(n1*n2),
                      gr13/sqrt(n1*n3),
                      gr14/sqrt(n1*n4),
                      gr12/sqrt(n1*n2),
                      1/n2,
                      gr23/sqrt(n2*n3),
                      gr24/sqrt(n2*n4),
                      gr13/sqrt(n1*n3),
                      gr23/sqrt(n2*n3),
                      1/n3,
                      gr34/sqrt(n3*n4),
                      gr14/sqrt(n1*n4),
                      gr24/sqrt(n2*n4),
                      gr34/sqrt(n3*n4),
                      1/n4)
  Sigma <- matrix(vecsigma,4,4)
  return(Sigma)
  
}
sigma = 0.4
gr12 = 0.8
gr13 = 0.8
gr14 = 0.8
gr23 = 0.8
gr24 = 0.8
gr34 = 0.8
n1 = n.EUR.snp
n2 = n.AFR.snp
n3 = n.AMR.snp
n4 = n.EAS.snp
# n1 = n.all.shared
# n2 = n.AFR.snp
# n3 = n.all.shared
# n4 = n.all.shared
Sigma <- GenSigma(sigma,n1,n2,n3,n4,
                  gr12,gr13,g14,
                  gr23,g24,g34)

save(Sigma, file ="/data/zhangh24/multi_ethnic/result/LD_simulation/causal_Sigma.rdata")
# beta_shared =  rmvnorm(n.all.shared,c(0,0,0,0),
#                 sigma=Sigma)
# beta[all.shared,] = beta_shared
# colnames(beta) <- paste0("beta_",c("EUR","AFR","AMR","EAS"))
# 
# AFR.specific = AFR.idx[AFR.idx%in%all.shared==F]
# n.AFR.specific = length(AFR.specific)
# beta_AFR_specific = rnorm(n.AFR.specific,mean = 0,
#                             sd = sqrt(sigma*1/n.AFR.snp))

#beta[AFR.specific,2] = beta_AFR_specific
#beta represent standarize scale effect-size
beta <- matrix(0,n.total.snp,4)
set.seed(666)
library(mvtnorm)
beta =  rmvnorm(n.total.snp,c(0,0,0,0),
                sigma=Sigma)
colnames(beta) <- paste0("beta_",c("EUR","AFR","AMR","EAS"))

cau.snp.infor <- cbind(cau.snp.infor,beta)
#transform the genfile to additive genotype file in EUR

n.snp <- length(which(cau.snp.infor$EUR>=0.01&
                        cau.snp.infor$EUR<=0.99))
n <- 120000

#gen.file <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/cau.combined.tag.gen"),header=F))
#snp.id <- as.character(gen.file[,3])
load("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/causal_genotype.rdata")
EUR.idx <- which(cau.snp.infor$EUR>=0.01&
                   cau.snp.infor$EUR<=0.99)

EUR.snp <- cau.snp.infor[EUR.idx,]

row.names(genotype) <- EUR.snp$id
beta_EUR <- EUR.snp$beta_EUR
MAF_EUR <- EUR.snp$EUR

sum(genotype[1,])/(2*ncol(genotype))
#MAF_test <- rowSums(genotype)/(2*ncol(genotype))


genotype_standard <- genotype

for(k in 1:nrow(genotype_standard)){
  genotype_standard[k,] = (genotype[k,]-2*MAF_EUR[k])/sqrt(2*MAF_EUR[k]*(1-MAF_EUR[k]))
}
G_value = t(genotype_standard)%*%(beta_EUR)
y = G_value+rnorm(n,mean=0,sd = sqrt(1-sigma))
mean(y)
var(y)

n.train = 100000
n.test <-20000
n.rep <- 100
y_test_mat <- matrix(0,n.test,n.rep)
for(l in 1:n.rep){
  y.temp = G_value+rnorm(n,mean=0,sd = sqrt(1-sigma))
  y_test_mat[,l] = y.temp[(n.train+1):length(y.temp)]
}
mean(y_test_mat[,1])
var(y_test_mat[,1])



#temp_var = apply(genotype_standard,1,var)
#beta_EUR_shared = beta_EUR[beta_EUR!=0]
#var(beta_EUR_shared)*length(beta_EUR_shared)

#sample <- read.table("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/sample.txt",header=T)


save(y,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/phenotype.rdata")
save(genotype,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/causal_genotype.rdata")
save(y_test_mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/phenotype_test_mat.rdata")





#transform the genfile to additive genotype file in AFR

n.snp <- length(which(cau.snp.infor$AFR>=0.01&
                        cau.snp.infor$AFR<=0.99))
load("/data/zhangh24/multi_ethnic/result/LD_simulation/AFR/causal_genotype.rdata")
AFR.idx <- which(cau.snp.infor$AFR>=0.01&
                   cau.snp.infor$AFR<=0.99)

AFR.snp <- cau.snp.infor[AFR.idx,]
# all.equal(paste0(AFR.snp$CHR,":",
#                  AFR.snp$position),
#           paste0(CHR,":",position))

row.names(genotype) <- AFR.snp$id
beta_AFR <- AFR.snp$beta_AFR
MAF_AFR <- AFR.snp$AFR

genotype_standard <- genotype
for(k in 1:nrow(genotype_standard)){
  genotype_standard[k,] = (genotype[k,]-2*MAF_AFR[k])/sqrt(2*MAF_AFR[k]*(1-MAF_AFR[k]))
}

n = ncol(genotype_standard)
G_value = t(genotype_standard)%*%(beta_AFR)
y = G_value+rnorm(n,mean=0,sd = sqrt(1-sigma))
mean(y)
var(y)

n.train = 15000
n.test <-3000
n.rep <- 100
y_test_mat <- matrix(0,n.test,n.rep)
for(l in 1:n.rep){
  y.temp = G_value+rnorm(n,mean=0,sd = sqrt(1-sigma))
  y_test_mat[,l] = y.temp[(n.train+1):length(y.temp)]
}
mean(y_test_mat[,2])
var(y_test_mat[,2])


save(y,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/AFR/phenotype.rdata")
save(genotype,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/AFR/causal_genotype.rdata")
save(y_test_mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/AFR/phenotype_test_mat.rdata")




#transform the genfile to additive genotype file in AMR

n.snp <- length(which(cau.snp.infor$AMR>=0.01&
                        cau.snp.infor$AMR<=0.99))
load("/data/zhangh24/multi_ethnic/result/LD_simulation/AMR/causal_genotype.rdata")
n <- 17830
AMR.idx <- which(cau.snp.infor$AMR>=0.01&
                   cau.snp.infor$AMR<=0.99)

AMR.snp <- cau.snp.infor[AMR.idx,]
beta_AMR <- AMR.snp$beta_AMR
MAF_AMR <- AMR.snp$AMR


genotype_standard <- genotype
for(k in 1:nrow(genotype_standard)){
  genotype_standard[k,] = (genotype[k,]-2*MAF_AMR[k])/sqrt(2*MAF_AMR[k]*(1-MAF_AMR[k]))
}
n = ncol(genotype_standard)
G_value = t(genotype_standard)%*%(beta_AMR)
y = G_value+rnorm(n,mean=0,sd = sqrt(1-sigma))
mean(y)
var(y)

n.train = 15000
n.test <-2830
n.rep <- 100
y_test_mat <- matrix(0,n.test,n.rep)
for(l in 1:n.rep){
  y.temp = G_value+rnorm(n,mean=0,sd = sqrt(1-sigma))
  y_test_mat[,l] = y.temp[(n.train+1):length(y.temp)]
}
mean(y_test_mat[,2])
var(y_test_mat[,2])



save(y,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/AMR/phenotype.rdata")
save(genotype,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/AMR/causal_genotype.rdata")
save(y_test_mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/AMR/phenotype_test_mat.rdata")



#transform the genfile to additive genotype file in EAS

n.snp <- length(which(cau.snp.infor$EAS>=0.01&
                        cau.snp.infor$EAS<=0.99))
load("/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/causal_genotype.rdata")
n <- 18000
EAS.idx <- which(cau.snp.infor$EAS>=0.01&
                   cau.snp.infor$EAS<=0.99)

EAS.snp <- cau.snp.infor[EAS.idx,]

beta_EAS <- EAS.snp$beta_EAS
MAF_EAS <- EAS.snp$EAS

genotype_standard <- genotype
for(k in 1:nrow(genotype_standard)){
  genotype_standard[k,] = (genotype[k,]-2*MAF_EAS[k])/sqrt(2*MAF_EAS[k]*(1-MAF_EAS[k]))
}

n = ncol(genotype_standard)
G_value = t(genotype_standard)%*%(beta_EAS)
y = G_value+rnorm(n,mean=0,sd = sqrt(1-sigma))
mean(y)
var(y)

n.train = 15000
n.test <-3000
n.rep <- 100
y_test_mat <- matrix(0,n.test,n.rep)
for(l in 1:n.rep){
  y.temp = G_value+rnorm(n,mean=0,sd = sqrt(1-sigma))
  y_test_mat[,l] = y.temp[(n.train+1):length(y.temp)]
}
mean(y_test_mat[,2])
var(y_test_mat[,2])

save(y,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/phenotype.rdata")
save(genotype,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/causal_genotype.rdata")
save(y_test_mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/phenotype_test_mat.rdata")


load("/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/phenotype.rdata")
load("/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/causal_genotype.rdata")

p_e <- rep(0,nrow(genotype))
for(k in 1:length(p_e)){
  print(k)
  model <- lm(y[1:n.train]~genotype[k,1:n.train])
  p_e[k] <- coefficients(summary(model))[2,4]
}

p_e_order <- p_e[order(p_e)]

p_order <- p[order(p)]
p_a_order <- p_a[order(p_a)]











sigma = 0.4
gr12 = 0.4
gr13 = 0.6
gr14 = 0.8
gr23 = 0.6
gr24 = 0.6
gr34 = 0.6
n1 = n.EUR.snp
n2 = n.AFR.snp
n3 = n.AMR.snp
n4 = n.EAS.snp

library(mvtnorm)
#beta represent standarize scale effect-size
set.seed(666)

#order the causal SNPs by CHR and position
cau.snp.infor <- cbind(cau.snp.infor,beta)
idx.order <- order(cau.snp.infor$CHR,
                   cau.snp.infor$position)
cau.snp.infor <- cau.snp.infor[idx.order,]
#save(cau.snp.infor,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/cau.snp.infor.effect.rdata")

#transform the genfile to additive genotype file in EUR and EAS

n.snp <- length(which(cau.snp.infor$EUR>=0.01&
                        cau.snp.infor$EUR<=0.99))
snp.id <- rep("c",n.snp)
position <- rep(0,n.snp)
CHR <- rep(0,n.snp)
n <- 120000
gen.impu <- matrix(0,n.snp,3*n)
total <- 0

#gen.file <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/cau.combined.tag.gen"),header=F))
#snp.id <- as.character(gen.file[,3])

for(j in 1:22){
  print(j)
  
  gen.file <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/chr",j,"_cau.combined.tag.gen"),header=F))
  temp <- nrow(gen.file)
  snp.id[total+(1:temp)] <- as.character(gen.file[,3])
  position[total+(1:temp)] <- as.numeric(gen.file[,4])
  CHR[total+(1:temp)] <- j
  gen.impu[total+(1:temp),] <- as.matrix(gen.file[,7:ncol(gen.file)])
  
  total <- total+temp
}
#calculate the genotype 0*P(aa)+1*P(Aa)+2*P(AA)
genotype <- matrix(0,n.snp,n)  
for(k in 1:n){
  if(k %%100==0){
    print(k)  
  }
  
  genotype[,k]<-  gen.impu[,3*k-1]+gen.impu[,3*k]*2
  
}
EUR.idx <- which(cau.snp.infor$EUR>=0.01&
                   cau.snp.infor$EUR<=0.99)

EUR.snp <- cau.snp.infor[EUR.idx,]
all.equal(paste0(EUR.snp$CHR,":",
                 EUR.snp$position),
          paste0(CHR,":",position))

row.names(genotype) <- EUR.snp$id
beta_EUR <- EUR.snp$beta_EUR
MAF_EUR <- EUR.snp$EUR

sum(genotype[1,])/(2*ncol(genotype))
#MAF_test <- rowSums(genotype)/(2*ncol(genotype))


genotype_standard <- genotype
for(k in 1:nrow(genotype_standard)){
  genotype_standard[k,] = (genotype[k,]-2*MAF_EUR[k])/sqrt(2*MAF_EUR[k]*(1-MAF_EUR[k]))
}
y = t(genotype_standard)%*%(beta_EUR)+rnorm(n,mean=0,sd = sqrt(1-sigma))
mean(y)
var(y)

sample <- read.table("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/sample.txt",header=T)


save(y,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/phenotype.rdata")
save(genotype,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/causal_genotype.rdata")






#transform the genfile to additive genotype file in AFR

n.snp <- length(which(cau.snp.infor$AFR>=0.01&
                        cau.snp.infor$AFR<=0.99))
snp.id <- rep("c",n.snp)
position <- rep(0,n.snp)
CHR <- rep(0,n.snp)
n <- 18000
gen.impu <- matrix(0,n.snp,3*n)
total <- 0

#gen.file <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/cau.combined.tag.gen"),header=F))
#snp.id <- as.character(gen.file[,3])

for(j in 1:22){
  print(j)
  
  gen.file <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/AFR/chr",j,"_cau.combined.tag.gen"),header=F))
  temp <- nrow(gen.file)
  snp.id[total+(1:temp)] <- as.character(gen.file[,3])
  position[total+(1:temp)] <- as.numeric(gen.file[,4])
  CHR[total+(1:temp)] <- j
  gen.impu[total+(1:temp),] <- as.matrix(gen.file[,7:ncol(gen.file)])
  
  total <- total+temp
}
#calculate the genotype 0*P(aa)+1*P(Aa)+2*P(AA)
genotype <- matrix(0,n.snp,n)  
for(k in 1:n){
  if(k %%100==0){
    print(k)  
  }
  
  genotype[,k]<-  gen.impu[,3*k-1]+gen.impu[,3*k]*2
  
}
AFR.idx <- which(cau.snp.infor$AFR>=0.01&
                   cau.snp.infor$AFR<=0.99)

AFR.snp <- cau.snp.infor[AFR.idx,]
all.equal(paste0(AFR.snp$CHR,":",
                 AFR.snp$position),
          paste0(CHR,":",position))

row.names(genotype) <- AFR.snp$id
beta_AFR <- AFR.snp$beta_AFR
MAF_AFR <- AFR.snp$AFR

genotype_standard <- genotype
for(k in 1:nrow(genotype_standard)){
  genotype_standard[k,] = (genotype[k,]-2*MAF_AFR[k])/sqrt(2*MAF_AFR[k]*(1-MAF_AFR[k]))
}


y = t(genotype_standard)%*%(beta_AFR)+rnorm(n,mean=0,sd = sqrt(1-sigma))
mean(y)
var(y)

save(y,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/AFR/phenotype.rdata")
save(genotype,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/AFR/causal_genotype.rdata")





#transform the genfile to additive genotype file in AMR

n.snp <- length(which(cau.snp.infor$AMR>=0.01&
                        cau.snp.infor$AMR<=0.99))
snp.id <- rep("c",n.snp)
position <- rep(0,n.snp)
CHR <- rep(0,n.snp)
n <- 17830
gen.impu <- matrix(0,n.snp,3*n)
total <- 0

#gen.file <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/cau.combined.tag.gen"),header=F))
#snp.id <- as.character(gen.file[,3])

for(j in 1:22){
  print(j)
  
  gen.file <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/AMR/chr",j,"_cau.combined.tag.gen"),header=F))
  temp <- nrow(gen.file)
  snp.id[total+(1:temp)] <- as.character(gen.file[,3])
  position[total+(1:temp)] <- as.numeric(gen.file[,4])
  CHR[total+(1:temp)] <- j
  gen.impu[total+(1:temp),] <- as.matrix(gen.file[,7:ncol(gen.file)])
  
  total <- total+temp
}
#calculate the genotype 0*P(aa)+1*P(Aa)+2*P(AA)
genotype <- matrix(0,n.snp,n)  
for(k in 1:n){
  if(k %%100==0){
    print(k)  
  }
  
  genotype[,k]<-  gen.impu[,3*k-1]+gen.impu[,3*k]*2
  
}
AMR.idx <- which(cau.snp.infor$AMR>=0.01&
                   cau.snp.infor$AMR<=0.99)

AMR.snp <- cau.snp.infor[AMR.idx,]
all.equal(paste0(AMR.snp$CHR,":",
                 AMR.snp$position),
          paste0(CHR,":",position))

row.names(genotype) <- AMR.snp$id
beta_AMR <- AMR.snp$beta_AMR
MAF_AMR <- AMR.snp$AMR


genotype_standard <- genotype
for(k in 1:nrow(genotype_standard)){
  genotype_standard[k,] = (genotype[k,]-2*MAF_AMR[k])/sqrt(2*MAF_AMR[k]*(1-MAF_AMR[k]))
}



y = t(genotype_standard)%*%(beta_AMR)+rnorm(n,mean=0,sd = sqrt(1-sigma))
mean(y)
var(y)
save(y,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/AMR/phenotype.rdata")
save(genotype,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/AMR/causal_genotype.rdata")




#transform the genfile to additive genotype file in EAS

n.snp <- length(which(cau.snp.infor$EAS>=0.01&
                        cau.snp.infor$EAS<=0.99))
snp.id <- rep("c",n.snp)
position <- rep(0,n.snp)
CHR <- rep(0,n.snp)
n <- 18000
gen.impu <- matrix(0,n.snp,3*n)
total <- 0

#gen.file <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/cau.combined.tag.gen"),header=F))
#snp.id <- as.character(gen.file[,3])

for(j in 1:22){
  print(j)
  
  gen.file <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/chr",j,"_cau.combined.tag.gen"),header=F))
  temp <- nrow(gen.file)
  snp.id[total+(1:temp)] <- as.character(gen.file[,3])
  position[total+(1:temp)] <- as.numeric(gen.file[,4])
  CHR[total+(1:temp)] <- j
  gen.impu[total+(1:temp),] <- as.matrix(gen.file[,7:ncol(gen.file)])
  
  total <- total+temp
}
#calculate the genotype 0*P(aa)+1*P(Aa)+2*P(AA)
genotype <- matrix(0,n.snp,n)  
for(k in 1:n){
  if(k %%100==0){
    print(k)  
  }
  
  genotype[,k]<-  gen.impu[,3*k-1]+gen.impu[,3*k]*2
  
}
EAS.idx <- which(cau.snp.infor$EAS>=0.01&
                   cau.snp.infor$EAS<=0.99)

EAS.snp <- cau.snp.infor[EAS.idx,]
all.equal(paste0(EAS.snp$CHR,":",
                 EAS.snp$position),
          paste0(CHR,":",position))

row.names(genotype) <- EAS.snp$id
beta_EAS <- EAS.snp$beta_EAS
MAF_EAS <- EAS.snp$EAS

genotype_standard <- genotype
for(k in 1:nrow(genotype_standard)){
  genotype_standard[k,] = (genotype[k,]-2*MAF_EAS[k])/sqrt(2*MAF_EAS[k]*(1-MAF_EAS[k]))
}


y = t(genotype_standard)%*%(beta_EAS)+rnorm(n,mean=0,sd = sqrt(1-sigma))
mean(y)
var(y)
save(y,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/phenotype.rdata")
save(genotype,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/causal_genotype.rdata")

load("/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/phenotype.rdata")
load("/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/causal_genotype.rdata")

p_e <- rep(0,nrow(genotype))
for(k in 1:length(p_e)){
  print(k)
  model <- lm(y[1:n.train]~genotype[k,1:n.train])
  p_e[k] <- coefficients(summary(model))[2,4]
}

p_e_order <- p_e[order(p_e)]

p_order <- p[order(p)]
p_a_order <- p_a[order(p_a)]