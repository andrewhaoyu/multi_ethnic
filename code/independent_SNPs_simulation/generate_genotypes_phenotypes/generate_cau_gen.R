#Goal: Generate cau genotyped data for simulation
#simulate phenotypes data for AFR, EUR, AMR
#MAF based on 1000KG
#sample size EUR n =120000
#sample size AFR n = 18000
#sample size AMR n = 18000
#heritability for EUR 0.5
#heritability for AFR 0.5
#heritability for AMR 0.5
#5000 causal SNPs for each population
#4000 shared causal SNPs
#Genetic correlation for the between EUR and AFR is 0.4
#Genetic correlation for the between EUR and LAC is 0.6
#Genetic correlation for the between LAC and AFR is 0.6
#1000 SNPs for each as independent causal
# load("/spin1/users/zhangh24/KG.vcf/MAF_result/pruned_MAF.Rdata")
#  set.seed(666)
#  colnames(pruned.snp.clean)[8] <- "MAF.LAC"
#  n.snp <- nrow(pruned.snp.clean)
#  pruned.snp.permu <- pruned.snp.clean[sample(c(1:n.snp),n.snp),]
#  save(pruned.snp.permu,file = "/spin1/users/zhangh24/KG.vcf/MAF_result/pruned_MAF_permu.Rdata")

load("/spin1/users/zhangh24/KG.vcf/MAF_result/pruned_MAF_permu.Rdata")
set.seed(666)
# idx <- which.min(pruned.snp.permu$MAF.AFR[1:5000])
# pruned.snp.permu[idx,]

n.snp <- nrow(pruned.snp.permu)
library(bc2)

GenearteGenotype <- function(n.sub,n.snp,MAF){
  result <- matrix(rbinom(n.sub*n.snp,2,MAF),n.sub,n.snp,byrow = T)
  return(result)
}
MAF.EUR <- pruned.snp.permu$MAF.EUR
MAF.AFR <- pruned.snp.permu$MAF.AFR
MAF.LAC <- pruned.snp.permu$MAF.LAC
n.cau <- 7000
n.EUR <- 120000
n.AFR <- 18000
n.AMR <- 18000
#since the SNP MAF data is permutated
#for EUR and AFR, we will use 1:4000 as shared SNPs
#use 4000:5000 as nonshared SNPs for EUR
#5000:6000 as nonshared SNPs for AFR
#7000:8000 as nonshared SNPs for AMR

genotype_EUR <- GenearteGenotype(n.EUR,n.cau, MAF.EUR[1:n.cau])
genotype_AFR <- GenearteGenotype(n.AFR,n.cau,MAF.AFR[1:n.cau])
genotype_LAC <- GenearteGenotype(n.AMR,n.cau,MAF.LAC[1:n.cau])

genotype <- list(genotype_EUR,
                 genotype_AFR,
                 genotype_LAC)
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis')
save(genotype,file = paste0("./multi_ethnic/result/pruned_geno/geno_",1))

#load(paste0("./multi_ethnic/result/pruned_geno/geno_",3))

#start.end <- startend(n.snp,23,1)












