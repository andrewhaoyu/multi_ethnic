

setwd("/data/zhangh24/test_dna")

#update the family id and individual id in the fam file
library(data.table)
fam = fread("chr22.fam")
for(i in 1:3){
  fam_sub = fam
  fam_sub[,c(1:2)] = as.integer((1:120000)+120000*(i-1))
  
  write.table(fam_sub, file = paste0("/data/zhangh24/test_dna/chr22_sub",i,".fam"),col.names = F, row.names = F, quote = F)
}
#some coding error with plink 1.9; it will write 100000 as 1e+05 in merge function
#to avoid the issue, we will just 
fam = data.frame(FID = as.character(c(1:300000)), IID = as.character(c(1:300000)))
write.table(fam, file = "/data/zhangh24/test_dna/keep_id",row.names = F, col.names = F, quote = F)
total <- 0
merge_list = rep("c",2)
temp = 1
for(i in 2:3){
  merge_list[temp] <- paste0("/data/zhangh24/test_dna/chr22_sub",i)
  temp = temp+1
}
write.table(merge_list,file = paste0("/data/zhangh24/test_dna/merge_list"),
            col.names = F,
            row.names = F,
            quote=F)
#merge all the seperate sample into one
#merge the seperate genotype data into one file
system(paste0("/data/zhangh24/software/plink2 ",
              "--bfile /data/zhangh24/test_dna/chr22_sub1 ",
              "--merge-list /data/zhangh24/test_dna/merge_list ",
              "--keep /data/zhangh24/test_dna/keep_id ",
              "--make-bed --out /data/zhangh24/test_dna/chr22_300k"))


pheno = data.frame(FID = c(1:300000), IID = c(1:300000), pheno1 = rbinom(300000,1,0.5)+1,pheno2 = rbinom(300000,1,0.5)+1,pheno3 = rbinom(300000,1,0.5)+1)
covar = data.frame(FID = c(1:300000), IID = c(1:300000), covar1 = rnorm(300000),
                   covar2 = rbinom(300000,1,0.5),covar3 = rbinom(300000,1,0.5),
                   covar4 = rbinom(300000,1,0.5),covar5 = rbinom(300000,1,0.5),
                   covar6 = rbinom(300000,1,0.5),covar7 = rbinom(300000,1,0.5),
                   covar8 = rbinom(300000,1,0.5),covar9 = rbinom(300000,1,0.5),
                   covar10 = rbinom(300000,1,0.5),covar11 = rbinom(300000,1,0.5))
write.table(pheno, file = "/data/zhangh24/test_dna/pheno_300k.phen", row.names = F, col.names = T, quote = F)
write.table(covar, file = "/data/zhangh24/test_dna/covar_300k.cov", row.names = F, col.names = T, quote = F)
