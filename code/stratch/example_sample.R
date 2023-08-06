#run GWAS using plink2
#snptest is very slow compared to plink2
#create the sample file matching the snptest requirement
#sample file contains ID_1,ID_2,missing,pheno1
#the first row represent variable status 0,0,0,P
n <- 1000000
#first create plain sample files for transforming gen data to bim data
ID_1 = ID_2 = paste0("AFR_", formatC(c(0,c(1:n)), format = "f", digits = 0, big.mark = ""))
missing = rep(0,n+1)
#only keep the training data
pheno1 = rep(NA,n+1)
pheno1[1] = "P"
pheno1[2:(n+1)] = rnorm(n)
pheno = data.frame(ID_1,ID_2,missing,pheno1)
head(pheno)
write.table(pheno,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/EUR/sample.txt",quote=F,col.names = T,row.names = F)
