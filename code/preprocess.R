#Goal: prepare the data from 1000 genome for the simulations
#data preprocess
#1. download the vcf data from 1000 genome
#2. use plink to do the LD pruning with 1000 genome 0.05 r2 and 1Mb

result <- matrix("c",22,1)
for(i in 1:22){
  result[i,1] <- paste0("/spin1/users/zhangh24/plink --vcf ALL.chr",i,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep /spin1/users/zhangh24/KG.vcf/ID.EUR --indep-pairwise 1000'kb' 1 0.1 --out /spin1/users/zhangh24/KG.vcf/prunded_result/chr_",i)
}
write.table(result,
  file="/spin1/users/zhangh24/KG.vcf/LD_pruned.sh",
  row.names = F,
  col.names = F,
  quote=F)
########read in phenotype information
library(data.table)
phenotype_ori <- as.data.frame(fread("/spin1/users/zhangh24/KG.vcf/integrated_call_samples_v2.20130502.ALL.ped"),stringasFactors=F)
phenotype_update <- as.data.frame(fread("/spin1/users/zhangh24/KG.vcf/integrated_call_samples_v3.20130502.ALL.panel"),stringasFactors=F)

phenotype <- merge(phenotype_update,phenotype_ori,
                  by.x="sample",
                  by.y = "Individual ID")
colnames(phenotype)[1] <- "Individual ID"



library(dplyr)
###########vcf transformed file will put individual id as family id, we can't find the family id from vcf files
CreateIDfile <- function(phenotype,ID){
  ID.file <- phenotype %>% 
    filter(super_pop==ID) %>% 
    select(c("Individual ID"))
  ID.result <- data.frame(ID.file,ID.file,stringsAsFactors = F)
  colnames(ID.result) <- c("Family ID",
                           "Individual ID")
  return(ID.result)
}

ID.EUR <- CreateIDfile(phenotype,"EUR")


ID.AFR <- CreateIDfile(phenotype,"AFR")

ID.AMR <- CreateIDfile(phenotype,"AMR")

write.table(ID.EUR,file = "/spin1/users/zhangh24/KG.vcf/ID.EUR",row.names = F,
            col.names = F,quote=F)
write.table(ID.AFR,file = "/spin1/users/zhangh24/KG.vcf/ID.AFR",row.names = F,
            col.names = F,quote=F)
write.table(ID.AMR,file = "/spin1/users/zhangh24/KG.vcf/ID.AMR",row.names = F,
            col.names = F,quote=F)


#####calculate the MAF for all of the SNPs for EUR and AFR
result <- matrix("c",66,1)
for(i in 1:22){
  result[i,1] <- paste0("/spin1/users/zhangh24/plink --vcf ALL.chr",i,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep /spin1/users/zhangh24/KG.vcf/ID.EUR --freq --out /spin1/users/zhangh24/KG.vcf/MAF_result/EUR_chr_",i)
}

for(i in 1:22){
  result[i+22,1] <- paste0("/spin1/users/zhangh24/plink --vcf ALL.chr",i,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep /spin1/users/zhangh24/KG.vcf/ID.AFR --freq --out /spin1/users/zhangh24/KG.vcf/MAF_result/AFR_chr_",i)
}
write.table(result,
            file="/spin1/users/zhangh24/KG.vcf/MAF_cal.sh",
            row.names = F,
            col.names = F,
            quote=F)

for(i in 1:22){
  result[i+44,1] <- paste0("/spin1/users/zhangh24/plink --vcf ALL.chr",i,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep /spin1/users/zhangh24/KG.vcf/ID.AMR --freq --out /spin1/users/zhangh24/KG.vcf/MAF_result/AMR_chr_",i)
}
write.table(result,
            file="/spin1/users/zhangh24/KG.vcf/MAF_cal.sh",
            row.names = F,
            col.names = F,
            quote=F)


#####merge all the maf files into one
##### head -1 chr_1.frq > all.freq 
##### tail -n +2 -q chr_*.frq >> all.freq
##### head -1 EUR_chr_1.frq >> all_EUR.freq
##### tail -n +2 -q EUR_chr_*.frq >> all_EUR.freq
##### head -1 AFR_chr_1.frq >> all_AFR.freq
##### tail -n +2 -q AFR_chr_*.frq >> all_AFR.freq
##### head -1 AMR_chr_1.frq >> all_AMR.freq
##### tail -n +2 -q AMR_chr_*.frq >> all_AMR.freq

#####read in the LD pruned SNPs
num.snp = rep(0,22)
for(i in 1:22){
  geno.file = paste0("/spin1/users/zhangh24/KG.vcf/prunded_result/chr_",i,".prune.in")
  temp = system(paste0("wc -l ",geno.file),intern=T)
  temp = as.numeric(gsub(geno.file,"",temp))
  num.snp[i] = temp
}
snp.pruned <- data.frame(rep("c",sum(num.snp)),stringsAsFactors = F)
library(data.table)
total <- 0
for(i in 1:22){
  snp.temp <- as.data.frame(fread(paste0("/spin1/users/zhangh24/KG.vcf/prunded_result/chr_",i,".prune.in"),header=F))
  snp.pruned[(total+1):(total+num.snp[i]),1] <- snp.temp
  total <- total+num.snp[i]
  
}


colnames(snp.pruned) <- "rs_id"

all.snp.EUR <- as.data.frame(fread("/spin1/users/zhangh24/KG.vcf/MAF_result/all_EUR.freq",header=T))
colnames(all.snp.EUR)[5] <- "MAF.EUR"

all.snp.AFR <-  as.data.frame(fread("/spin1/users/zhangh24/KG.vcf/MAF_result/all_AFR.freq",header=T))
all.snp.AMR <- as.data.frame(fread("/spin1/users/zhangh24/KG.vcf/MAF_result/all_AMR.freq",header=T))
colnames(all.snp.AFR)[5] <- "MAF.AFR"
colnames(all.snp.AMR)[5] <- "MAF.AMR"
all.equal(all.snp.AFR$SNP,all.snp.EUR$SNP)
all.equal(all.snp.AFR$SNP,all.snp.AMR$SNP)

all.snp <- data.frame(all.snp.EUR,all.snp.AFR$MAF.AFR,all.snp.AMR$MAF.AMR)

dim(all.snp)
pruned.snp.infor <- merge(snp.pruned,all.snp,
                          by.x = "rs_id",
                          by.y = "SNP")
colnames(pruned.snp.infor)[c(7,8)] <- c("MAF.AFR","MAF.AMR")

#library(dplyr)
pruned.snp.clean= pruned.snp.infor %>% 
  filter(MAF.EUR>=0.05&
        rs_id!="."&
        MAF.AFR!=0&
        MAF.AMR!=0)
save(pruned.snp.clean,file= "/spin1/users/zhangh24/KG.vcf/MAF_result/pruned_MAF.Rdata")
colnames(all.snp)[c(7,8)] <- c("MAF.AFR",
                               "MAF.AMR")
all.snp <- all.snp[,-6]
all.snp.update = all.snp %>% 
  filter(MAF.EUR>=0.05)
save(all.snp.update,file = "/spin1/users/zhangh24/KG.vcf/MAF_result/all_snp_MAF.Rdata")
# all.snp.update = %>% all.snp %>% 
#   select(c(CHR,SNP,A1,A2,MAF.EUR,
#            all.snp.AFR.MAF.AFR))
#load("/spin1/users/zhangh24/KG.vcf/MAF_result/pruned_MAF.Rdata")

# x= matrix(rnorm(1000),ncol=1)
# x = rep(1,1000)
# y =x+ rnorm(1000)
# 
# x2 <- cbind(1,1,x)
# 
# 
# library(microbenchmark)
# ?microbenchmark
# 
# microbenchmark(fastLm(X=cbind(1,x),y=y),lm(y~x),
#                fastLm(X=x2[,c(1,3)],y=y))
# 
# model = fastLm(X=cbind(1,x),y=y)
# model1 <- fastLm(X=cbind(1,x),y=y)
# model2 <- lm(y~x)
# coef((model))
# install.packages("RcppArmadillo")
# library(RcppArmadillo)
# ?fastLm()


