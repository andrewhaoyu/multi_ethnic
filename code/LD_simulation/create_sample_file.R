#run GWAS using plink2
#snptest is very slow compared to plink2
#create the sample file matching the snptest requirement
#load the EUR phenotypes file
load("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/phenotype.rdata")
#phenotype format follows snptest format
#sample file contains ID_1,ID_2,missing,pheno1
#the first row represent variable status 0,0,0,P
n <- length(y)
ID_1 = c(0,c(1:n))
ID_2  = c(0,c(1:n))
missing = rep(0,n+1)
#only keep the training data
n.train = 100000
pheno1 = rep(NA,n+1)
pheno1[1] = "P"
pheno1[2:(n.train+1)] = y[1:n.train]
pheno = data.frame(ID_1,ID_2,missing,pheno1)
head(pheno)
write.table(pheno,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/sample.txt",quote=F,col.names = T,row.names = F)


#create the phenotypes file matching the snptest requirement
#load the AFR phenotypes file
load("/data/zhangh24/multi_ethnic/result/LD_simulation/AFR/phenotype.rdata")
#phenotype format follows snptest format
#sample file contains ID_1,ID_2,missing,pheno1
#the first row represent variable status 0,0,0,P
n <- length(y)
ID_1 = c(0,c(1:n))
ID_2  = c(0,c(1:n))
missing = rep(0,n+1)
#only keep the training data
n.train = 15000
pheno1 = rep(NA,n+1)
pheno1[1] = "P"
pheno1[2:(n.train+1)] = y[1:n.train]
pheno = data.frame(ID_1,ID_2,missing,pheno1)
head(pheno)
write.table(pheno,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/AFR/sample.txt",quote=F,col.names = T,row.names = F)



#create the phenotypes file matching the snptest requirement
#load the AMR phenotypes file
load("/data/zhangh24/multi_ethnic/result/LD_simulation/AMR/phenotype.rdata")
#phenotype format follows snptest format
#sample file contains ID_1,ID_2,missing,pheno1
#the first row represent variable status 0,0,0,P
n <- length(y)
ID_1 = c(0,c(1:n))
ID_2  = c(0,c(1:n))
missing = rep(0,n+1)
#only keep the training data
n.train = 15000
pheno1 = rep(NA,n+1)
pheno1[1] = "P"
pheno1[2:(n.train+1)] = y[1:n.train]
pheno = data.frame(ID_1,ID_2,missing,pheno1)
head(pheno)
write.table(pheno,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/AMR/sample.txt",quote=F,col.names = T,row.names = F)


#create the phenotypes file matching the snptest requirement
#load the EAS phenotypes file
load("/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/phenotype.rdata")
#phenotype format follows snptest format
#sample file contains ID_1,ID_2,missing,pheno1
#the first row represent variable status 0,0,0,P
n <- length(y)
ID_1 = c(0,c(1:n))
ID_2  = c(0,c(1:n))
missing = rep(0,n+1)
#only keep the training data
n.train = 15000
pheno1 = rep(NA,n+1)
pheno1[1] = "P"
pheno1[2:(n.train+1)] = y[1:n.train]
pheno = data.frame(ID_1,ID_2,missing,pheno1)
head(pheno)
write.table(pheno,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/sample.txt",quote=F,col.names = T,row.names = F)







library(data.table)
library(dplyr)
sample.temp <- as.data.frame(fread("/data/zhangh24/software/snptest/example/cohort1.sample",header=T))

sample.new = sample.temp %>% select(ID_1,ID_2,missing,pheno1)
sample.new$pheno1[2:nrow(sample.new)] <- rnorm(nrow(sample.new)-1,0,1)
write.table(sample.new,file = "/data/zhangh24/software/snptest/example/cohort1_test.sample",col.names = T,row.names = F,quote=F)

/data/zhangh24/software/snptest/snptest -data /data/zhangh24/multi_ethnic/result/LD_simulation/EUR/chr22.combined.tag.gen /data/zhangh24/multi_ethnic/result/LD_simulation/EUR/sample.txt -o /data/zhangh24/multi_ethnic/result/LD_simulation/EUR/summary.out -frequentist 1 -method score -pheno pheno1

/data/zhangh24/software/temp/plink2 --gen /data/zhangh24/multi_ethnic/result/LD_simulation/EUR/test_new.tag.gen --sample /data/zhangh24/multi_ethnic/result/LD_simulation/EUR/sample.txt --out /data/zhangh24/multi_ethnic/result/LD_simulation/EUR/summary_temp.out --linear 




plink.result <- as.data.frame(fread("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/summary_temp.out.pheno1.glm.linear",header=T))


load("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/causal_genotype.rdata")
load("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/phenotype.rdata")


model <- lm(y[1:n.train]~genotype[5000,1:n.train])
summary(model)

idx <- which(plink.result$ID=="rs7410305:50850541:G:A")
plink.result[idx,]
