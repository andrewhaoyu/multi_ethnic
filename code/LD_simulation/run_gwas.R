#run GWAS using plink2
#snptest is very slow compared to plink2
#create the phenotypes file matching the snptest requirement
#use plink2 to run


code <- rep("c",1000)
temp <- 1
eth <- c("EUR","AFR","AMR","EAS")
for(i in 1:4){
  for(j in 1:22){
    code[temp] <- paste0("/data/zhangh24/software/plink2 --gen /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/chr",j,".plink.tag.gen --sample /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/sample.txt --out /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary_chr",j,".out --linear ")
    temp <- temp+1
  }
}

code <- code[1:(temp-1)]
write.table(code,file = "/data/zhangh24/multi_ethnic/code/LD_simulation/run_gwas.sh",quote=F,row.names = F,col.names = F)











# 
# /data/zhangh24/software/snptest/snptest -data /data/zhangh24/multi_ethnic/result/LD_simulation/EUR/chr22.combined.tag.gen /data/zhangh24/multi_ethnic/result/LD_simulation/EUR/sample.txt -o /data/zhangh24/multi_ethnic/result/LD_simulation/EUR/summary.out -frequentist 1 -method score -pheno pheno1
i <-2
eth <- c("EUR","AFR","AMR","EAS")
pheno <- fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/sample.txt"))

idx <- which(is.na(pheno$pheno1))


n.train.vec <- c(100000,15000,15000,15000)


plink.result <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out"),header=T))





load(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/causal_genotype.rdata"))
load(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/phenotype.rdata"))

n.train <- n.train.vec[i]
model <- lm(y[1:n.train]~genotype[2500,1:n.train])
summary(model)
row.names(genotype)[2500]
idx <- which(plink.result$ID==row.names(genotype)[2500])
plink.result[idx,]







test.genotype <- genotype[,c(180*(c(1:100)-1)+1)]
cor.temp <- cor(test.genotype)
for(i in 1:100){
  for(j in 1:100){
    if(i>j&cor.temp[i,j]>=0.85)
      print(c(i,j))
  }
}

library(data.table)
#load the example file
example.gen <- as.data.frame(fread("/data/zhangh24/software/temp/example/ex.out.controls.gen",header=F))
genotype <- matrix(0,1000,100)
#genotype matrix with additive genotype value
#row represent SNP
#column represent subject
for(k in 1:100){
  genotype[,k] <- example.gen[,3*(k-1)+7]+2*example.gen[,3*(k-1)+8]
}
#the first 5 subjects correlation 
cor(genotype[,1:5])
