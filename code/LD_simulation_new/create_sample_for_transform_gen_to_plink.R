#run GWAS using plink2
#snptest is very slow compared to plink2
#create the sample file matching the snptest requirement
#load the EUR phenotypes file
#load("/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/phenotype.rdata")
#phenotype format follows snptest format
#sample file contains ID_1,ID_2,missing,pheno1
#the first row represent variable status 0,0,0,P
#n <- length(y)
#first create plain sample files for transforming gen data to bim data
n <- 120
ID_1 = c(0,c(1:n))
ID_2  = c(0,c(1:n))
missing = rep(0,n+1)
#only keep the training data
pheno1 = rep(NA,n+1)
pheno1[1] = "P"
pheno1[2:n] = NA
pheno = data.frame(ID_1,ID_2,missing,pheno1)
head(pheno)
write.table(pheno,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/sample_small.txt",quote=F,col.names = T,row.names = F)


#create the phenotypes file matching the snptest requirement
#load the AFR phenotypes file
#load("/data/zhangh24/multi_ethnic/result/LD_simulation/AFR/phenotype.rdata")
#phenotype format follows snptest format
#sample file contains ID_1,ID_2,missing,pheno1
#the first row represent variable status 0,0,0,P
n <- 120
ID_1 = c(0,c(1:n))
ID_2  = c(0,c(1:n))
missing = rep(0,n+1)
#only keep the training data
pheno1 = rep(NA,n+1)
pheno1[1] = "P"
pheno = data.frame(ID_1,ID_2,missing,pheno1)
head(pheno)
write.table(pheno,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/AFR/sample_small.txt",quote=F,col.names = T,row.names = F)



#create the phenotypes file matching the snptest requirement
#load the AMR phenotypes file
#phenotype format follows snptest format
#sample file contains ID_1,ID_2,missing,pheno1
#the first row represent variable status 0,0,0,P
n <- 120
ID_1 = c(0,c(1:n))
ID_2  = c(0,c(1:n))
missing = rep(0,n+1)
#only keep the training data
pheno1 = rep(NA,n+1)
pheno1[1] = "P"
pheno = data.frame(ID_1,ID_2,missing,pheno1)
head(pheno)
write.table(pheno,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/AMR/sample_small.txt",quote=F,col.names = T,row.names = F)


#create the phenotypes file matching the snptest requirement
#load the EAS phenotypes file
#phenotype format follows snptest format
#sample file contains ID_1,ID_2,missing,pheno1
#the first row represent variable status 0,0,0,P
n <- 120
ID_1 = c(0,c(1:n))
ID_2  = c(0,c(1:n))
missing = rep(0,n+1)
#only keep the training data
pheno1 = rep(NA,n+1)
pheno1[1] = "P"
pheno = data.frame(ID_1,ID_2,missing,pheno1)
head(pheno)
write.table(pheno,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EAS/sample_small.txt",quote=F,col.names = T,row.names = F)



n <- 120
ID_1 = c(0,c(1:n))
ID_2  = c(0,c(1:n))
missing = rep(0,n+1)
#only keep the training data
pheno1 = rep(NA,n+1)
pheno1[1] = "P"
pheno = data.frame(ID_1,ID_2,missing,pheno1)
head(pheno)
write.table(pheno,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/SAS/sample_small.txt",quote=F,col.names = T,row.names = F)
