#Goal: use R to seperate haplotype data by ethnic groups
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
i2 = as.numeric(args[[2]])
library(dplyr)
library(data.table)
sample <- read.table("/spin1/users/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3.sample",header=T)
country <- names(table(sample$GROUP))

#hap file is coded with row as SNP, columns as subjects
#two continious columns represent one subject
library(data.table)
haplotype <- as.data.frame(fread("/spin1/users/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",i2,".hap",header=F))

#find all the id in specific country group
idx <- which(sample$GROUP%in%country[i1])

idx.sub <- rep(0,length(idx))
temp <- 1
for(i in 1:length(idx)){
  idx.sub[temp] <- 2*idx[i]-1
  temp <- temp+1
  idx.sub[temp] <- 2*idx[i]
  temp <- temp+1
}

hap.sub <- haplotype[,idx.sub,drop=F]
write.table(hap.sub,file =paste0("/spin1/users/zhangh24/KG.impute2/",country[i1],"/chr",i2,".hap"),
            row.names = F,
            col.names=F,
            quote=F)