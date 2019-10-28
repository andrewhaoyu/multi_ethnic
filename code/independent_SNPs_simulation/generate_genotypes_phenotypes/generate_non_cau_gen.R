#Goal: Generate non causal genotyped data for simulation
#simulate phenotypes data for AFR, EUR, LAC
#MAF based on 1000KG
#sample size EUR n =120000
#sample size AFR n = 18000
#sample size LAC n = 18000
#heritability for EUR 0.5
#heritability for AFR 0.5
#heritability for LAC 0.5ls
#5000 causal SNPs for each population
#4000 shared causal SNPs
#Genetic correlation for the between EUR and AFR is 0.4
#Genetic correlation for the between EUR and LAC is 0.6
#Genetic correlation for the between LAC and AFR is 0.6
#1000 SNPs for each as independent causal
# load("/spin1/users/zhangh24/KG.vcf/MAF_result/pruned_MAF.Rdata")
#  set.seed(666)
#  n.snp <- nrow(pruned.snp.clean)
#  pruned.snp.permu <- pruned.snp.clean[sample(c(1:n.snp),n.snp),]
#  save(pruned.snp.permu,file = "/spin1/users/zhangh24/KG.vcf/MAF_result/pruned_MAF_permu.Rdata")
arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])

load("/spin1/users/zhangh24/KG.vcf/MAF_result/pruned_MAF_permu.Rdata")
set.seed(i1)
# idx <- which.min(pruned.snp.permu$MAF.AFR[1:5000])
# pruned.snp.permu[idx,]

n.snp <- nrow(pruned.snp.permu)
library(bc2)

GenearteGenotype <- function(n.sub,n.snp,MAF){
  result <- matrix(rbinom(n.sub*n.snp,2,MAF),n.sub,n.snp,byrow = T)
  return(result)
}
#split the genotype data into 500 different files
#1 is the causal genotype
#2:499 are the noncausal genotype
size <- 499
MAF.EUR <- pruned.snp.permu$MAF.EUR
MAF.AFR <- pruned.snp.permu$MAF.AFR
MAF.LAC <- pruned.snp.permu$MAF.LAC
n.cau <- 7000
start.end <- startend(n.snp-n.cau,size,i1)
start <- start.end[1]+n.cau
end <- start.end[2]+n.cau
n.EUR <- 120000
n.AFR <- 18000
n.LAC <- 18000
#since the SNP MAF data is permutated
#for EUR and AFR, we will use 1:4000 as shared SNPs
#use 4000:5000 as nonshared SNPs for EUR
#5000:6000 as nonshared SNPs for AFR
#6000:7000 as nonshared SNPs for LAC
n.noncau <- end-start+1
genotype_EUR <- GenearteGenotype(n.EUR,n.noncau, MAF.EUR[start:end])
genotype_AFR <- GenearteGenotype(n.AFR,n.noncau,MAF.AFR[start:end])
genotype_LAC <- GenearteGenotype(n.LAC,n.noncau,MAF.LAC[start:end])

genotype <- list(genotype_EUR,
                 genotype_AFR,
                 genotype_LAC)
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis')
save(genotype,file = paste0("./multi_ethnic/result/pruned_geno/geno_",i1+1))



