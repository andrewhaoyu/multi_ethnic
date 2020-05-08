library(data.table)
setwd("/Users/haoyuzhang/GoogleDrive/multi_ethnic")

#ex.haps <- as.data.frame(read.table(gzfile("/dcl01/chatterj/data/hzhang1/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz")))

ex.haps <- as.data.frame(read.table("./Hapgen2_test/1kg/test.hap",header=F))

n.sub <- ncol(ex.haps)/2
n.snp <- nrow(ex.haps)

genotype <- matrix(0,n.snp,n.sub)

for(k in 1:n.sub){
  genotype[,k] <- ex.haps[,2*k-1]+ex.haps[,2*k]
}

MAF <- rowSums(genotype)/(2*n.sub)

idx <- which(MAF>=0.05&MAF<=0.95)
genotype_common <- genotype[idx,]

genotype_new <- genotype_common
for(k in 1:length(idx)){
  genotype_new[k,] <- (genotype_common[k,]-2*MAF[idx[k]])/(sqrt(2*MAF[idx[k]]*(1-MAF[idx[k]])))
}

library(ggplot2)

data <- data.frame((cor(genotype_new)[upper.tri(cor(genotype_new))]))
colnames(data) <- "correlation"

ggplot(data,aes(correlation))+geom_histogram()+
  theme_Publication()+
  xlab("genotypic correlation")
range(data)









ex.gen <- as.data.frame(fread("./Hapgen2_test/hapgen2/example/ex.out.controls.gen"))


n.sub <- (ncol(ex.gen)-5)/3
n.snp <- nrow(ex.gen)
genotype <- matrix(0,n.snp,n.sub)
for(k in 1:n.sub){
  genotype[,k] <- ex.gen[,3*k+4]+2*ex.gen[,3*k+5]
}
cor(genotype[,1],genotype[,2])
MAF <- rowSums(genotype)/(2*n.sub)
idx <- which(MAF>=0.05&MAF<=0.95)
genotype_common <- genotype[idx,]
genotype_new <- genotype_common
for(k in 1:length(idx)){
  genotype_new[k,] <- (genotype_common[k,]-2*MAF[idx[k]])/(sqrt(2*MAF[idx[k]]*(1-MAF[idx[k]])))
}


n.sub <- ncol(ex.haps)/2
n.snp <- nrow(ex.haps)

library(ggplot2)

data <- data.frame((cor(genotype_new[,1:20])[upper.tri(cor(genotype_new[,1:20]))]))
colnames(data) <- "correlation"
library(ggthems)
ggplot(data,aes(correlation))+geom_histogram()+
  theme_Publication()+
  xlab("genotypic correlation")
range(data)

ex.haps <- as.data.frame(fread("./Hapgen2_test/hapgen2/example/ex.haps"))
ex.hap1 <- ex.haps[,1:60]
MAF <- rowSums(ex.hap1)/(2*ncol(ex.hap1))

ex.hap2 <- ex.haps[,61:120]
write.table(ex.hap1,file = paste0("./Hapgen2_test/hapgen2/example/ex.hap1"),row.names = F,col.names = F,quote = F)
write.table(ex.hap2,file = paste0("./Hapgen2_test/hapgen2/example/ex.hap2"),row.names = F,col.names = F,quote = F)

ex.leg <- as.data.frame(fread("./Hapgen2_test/hapgen2/example/ex.leg",header = T))





ex.gen <- as.data.frame(fread("./Hapgen2_test/hapgen2/example/ex.out1.controls.gen"))


n.sub <- (ncol(ex.gen)-5)/3
n.snp <- nrow(ex.gen)
genotype <- matrix(0,n.snp,n.sub)
for(k in 1:n.sub){
  genotype[,k] <- ex.gen[,3*k+4]+2*ex.gen[,3*k+5]
}

MAF <- rowSums(genotype)/(2*n.sub)
idx <- which(MAF>=0.05&MAF<=0.95)
genotype_common <- genotype[idx,]
genotype_new <- genotype_common
for(k in 1:length(idx)){
  genotype_new[k,] <- (genotype_common[k,]-2*MAF[idx[k]])/(sqrt(2*MAF[idx[k]]*(1-MAF[idx[k]])))
}

genotype_new1 = genotype_new


ex.gen <- as.data.frame(fread("./Hapgen2_test/hapgen2/example/ex.out2.controls.gen"))


n.sub <- (ncol(ex.gen)-5)/3
n.snp <- nrow(ex.gen)
genotype <- matrix(0,n.snp,n.sub)
for(k in 1:n.sub){
  genotype[,k] <- ex.gen[,3*k+4]+2*ex.gen[,3*k+5]
}

# MAF <- rowSums(genotype)/(2*n.sub)
# idx <- which(MAF>=0.05&MAF<=0.95)
genotype_common <- genotype[idx,]
genotype_new <- genotype_common
for(k in 1:length(idx)){
  genotype_new[k,] <- (genotype_common[k,]-2*MAF[idx[k]])/(sqrt(2*MAF[idx[k]]*(1-MAF[idx[k]])))
}

genotype_new2 = genotype_new
dim(genotype_new1)
dim(genotype_new2)
#range(cor(genotype_new1)[upper.tri(cor(genotype_new1))])
dim(cor(genotype_new1,genotype_new2))
range((cor(genotype_new1,genotype_new2)))

data <- data.frame(as.vector(cor(genotype_new1,genotype_new2)))
colnames(data) <- "correlation"

ggplot(data,aes(correlation))+geom_histogram()+
  theme_Publication()+
  xlab("genotypic correlation")
range(data)



n.snp <- 487
n.sub <- 100
genotype <- matrix(0,n.snp,n.sub)
genotype_new <- matrix(0,n.snp,n.sub)
MAF <- runif(487,0.05,0.95)
for(k in 1:n.snp){
  genotype[k,] <- rbinom(n.sub,2,MAF[k])
  genotype_new[k,] <- (genotype[k,]-MAF[k]*2)/(sqrt(2*MAF[k]*(1-MAF[k])))
}
range(cor(genotype_new)[upper.tri(cor(genotype_new))])




    
ex.haps <- as.data.frame(fread("./Hapgen2_test/1kg/test.hap"),header=F)  
dim(ex.haps)
n.sub <- ncol(ex.haps)/2
n.snp <- nrow(ex.haps)

genotype <- matrix(0,n.snp,n.sub)

for(k in 1:n.sub){
  genotype[,k] <- ex.haps[,2*k-1]+ex.haps[,2*k]
}

sample <- as.data.frame(fread("./Hapgen2_test/1kg/1000GP_Phase3.sample",header=T))

MAF <- rowSums(genotype)/(2*n.sub)

idx <- which(MAF>=0.05&MAF<=0.95)
genotype_common <- genotype[idx,]
length(idx)
genotype_new <- genotype_common
for(k in 1:length(idx)){
  genotype_new[k,] <- (genotype_common[k,]-2*MAF[idx[k]])/(sqrt(2*MAF[idx[k]]*(1-MAF[idx[k]])))
}
jdx <- which(sample$POP=="GBR")
genotype_new_tar <- genotype_new[,jdx]

idx1 <- which(sample$ID=="HG03343"
                )
idx2 <- which(sample$ID=="HG03352")
sample[idx1,]
sample[idx2,]
cor(genotype_common[,idx1],
    genotype_common[,idx2])

data <- data.frame((cor(genotype_new_tar)[upper.tri(cor(genotype_new_tar))]))
colnames(data) <- "correlation"
range(data)
which(cor(genotype_new_tar)<=0.53&
               cor(genotype_new_tar)>=0.52,arr.ind = T)
sample[23,]
sample[19,]

range(cor(genotype_new[,1:100])[upper.tri(cor(genotype_new[,1:100]))])
cor(genotype_new[,1:100])[upper.tri(cor(genotype_new[,1:100]))][250]
a = cor(genotype_new[,1:100])
which(0.5416110<a&a<0.5416112,arr.ind = T)
ggplot(data,aes(correlation))+geom_histogram()+
  theme_Publication()+
  xlab("genotypic correlation")
range(data)


sample <- as.data.frame(fread("./Hapgen2_test/1kg/1000GP_Phase3.sample",header=T))
idx <- which(sample$ID=="HG00405")
idx
