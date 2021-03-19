#LD_clumping for ARIC data
#load the p-value results
args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

#i = as.numeric(args[[1]])
#l = as.numeric(args[[2]])
i = 1

l = 1
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)

eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")
library(bigsnpr)
# for(i in 1:2){
#   for(l in c(1,3)){
setwd("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic")
temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")


pheno <- as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/pheno/pheno.txt")))
colnames(pheno)[2] = "ID"
genotype.fam <- as.data.frame(fread(paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc1.fam")))
load(paste0(temp.dir,"multi_PRS.rdata"))
prs.data = data.frame(ID = genotype.fam$V2,prs= multi_PRS,stringsAsFactors = F)
prs.all = left_join(pheno,prs.data,by="ID")



startend <- function(num,size,ind){
  split.all <- split(1:num,cut(1:num,size))
  temp <- split.all[[ind]]
  start <- temp[1]
  end <- temp[length(temp)]
  return(c(start,end))
}
n.rep = 5
rer2_prs_pc_rep =rep(0,n.rep)
#split the data into 5 fold
for(i_rep in 1:n.rep){
  start.end <- startend(nrow(prs.all),n.rep,i_rep)
  vad.id = c(start.end[1]:start.end[2])
  test.id = setdiff(c(1:nrow(prs.all)),vad.id)
  
  prs.score = prs.all[,16:ncol(prs.all)]
  prs.score = as_FBM(prs.score)
  model1.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.all)
  y = model1.null$residual
  penalized_model = big_spLinReg(prs.score,y,alphas = c(1, 0.01, 1e-04),K=4)
  
  
  
  
  
  prs.test = prs.all[test.id,]
  prs.vad = prs.all[vad.id,]

  prs.score.test = prs.test[,16:ncol(prs.test)]
  prs.score.test = as_FBM(prs.score.test)
  model1.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.test)
  y = model1.null$residual
  penalized_model = big_spLinReg(prs.score.test,y,alphas = c(1, 0.01, 1e-04),K=4)
  idx <- which.min(summary(penalized_model)$validation_loss)
  new_beta = summary(penalized_model)$beta[[idx]]
  
  prs.score.vad = as.matrix(prs.vad[,16:ncol(prs.test)])
  prs.score.vad.com <- prs.score.vad%*%new_beta
  model2.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.vad)
  model2.prs <- lm(model2.null$residual~prs.score.vad.com)
  rer2_prs_pc_rep[i_rep] = summary(model2.prs)$r.square
}

# multi_PRS_new = as_FBM(multi_PRS_new)
# y = rnorm(9345)
# temp = big_spLinReg(multi_PRS_new,y,alphas = c(1, 0.01, 1e-04))
