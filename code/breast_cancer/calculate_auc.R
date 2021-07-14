startend <- function(num,size,ind){
  split.all <- split(1:num,cut(1:num,size))
  temp <- split.all[[ind]]
  start <- temp[1]
  end <- temp[length(temp)]
  return(c(start,end))
}

args = commandArgs(trailingOnly = T)
l = as.numeric(args[[1]])
i1 = as.numeric(args[[2]])
library(tidyverse)
library(data.table)
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5,1.0)
setwd("/data/zhangh24/multi_ethnic/data/")
library(data.table)
library(tidyverse)
library(RISCA)
library(boot)
trait = c("overall","erpos","erneg")
data.dir = "/data/zhangh24/multi_ethnic/data/"
l = 1
out.dir.prs = "/data/zhangh24/multi_ethnic/result/breast_cancer/prs/"
n.rep = 2

auc.vec.test.rep = matrix(0,length(pthres),n.rep)
auc.vec.vad.rep = matrix(0,length(pthres),n.rep)

for(k in 1:length(pthres)){
  if(i1==1){
    prs.score = fread(paste0(out.dir.prs,"prs_",trait[l],".p_value_",k,".profile"),header=T)
  }else{
    prs.score = fread(paste0(out.dir.prs,"prs_",trait[l],"_mega.p_value_",k,".profile"),header=T)
  }
  
  
  
  prs.score = prs.score %>% 
    separate(IID,into=c("ID","ohter_id","nci_id"),remove=F)
  pheno <- fread(paste0(data.dir,"GBHS_pheno.csv"))
  colnames(pheno)[5] = "ER_status"
  
  pheno = pheno %>% 
    mutate(overall_status=ifelse(Status=="Control",0,1))
  
  pheno.update = left_join(pheno,prs.score,by=c("SB_ID"="ID"))
  #cross-validation for the data
  for(i_rep in 1:n.rep){
    start.end <- startend(nrow(pheno.update),n.rep,i_rep)
    vad.id = c(start.end[1]:start.end[2])
    test.id = setdiff(c(1:nrow(pheno.update)),vad.id)
    test.data <- pheno.update[test.id,]
    vad.data <- pheno.update[vad.id,]
    roc_obj <- roc.binary(status="overall_status", variable=paste0("SCORE",k,"_AVG"),
                          confounders=~EV1+EV2+EV3+EV4+EV5+
                            EV6+EV7+EV8+EV9+EV10,
                          data=test.data, precision=seq(0.1,0.9, by=0.1))
    auc.vec.test.rep[k,i_rep] = roc_obj$auc
    roc_obj <- roc.binary(status="overall_status", variable=paste0("SCORE",k,"_AVG"),
                          confounders=~EV1+EV2+EV3+EV4+EV5+
                            EV6+EV7+EV8+EV9+EV10,
                          data=vad.data, precision=seq(0.1,0.9, by=0.1))
    
    auc.vec.vad.rep[k,i_rep] = roc_obj$auc
  }
  
  }
  
  
  
clump.snp <- LD[,3,drop=F] 
sum.data = sum.data %>% 
  select(CHR,ID,POS,Effect_allele,Effect,P) %>% 
  rename(SNP=ID,
         A1 = Effect_allele,
         BETA = Effect)

prs.all <- left_join(clump.snp,sum.data,by="SNP")
#find all the duplicated SNPs and keep the more significant one
dup.id <- prs.all$SNP[duplicated(prs.all$SNP)]
if(length(dup.id)!=0){
  remove.idx = rep(0,length(dup.id))
  for(k in 1:length(dup.id)){
    jdx <- which(prs.all$SNP==dup.id[k])
    which.max(prs.all$P[jdx])
    remove.idx[k] = jdx[which.max(prs.all$P[jdx])]
  }
  prs.all = prs.all[-remove.idx,]
  
}
n_pthres <- length(pthres)
q_range = data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres))
write.table(q_range,file = paste0(temp.dir,"q_range_file"),col.names = T,row.names = F,quote=F)
prs.file <- prs.all %>% 
  select(SNP,A1,BETA)
write.table(prs.file,file = paste0(temp.dir,"prs_file"),col.names = T,row.names = F,quote=F)
p.value.file <- prs.all %>% 
  select(SNP,P)
write.table(p.value.file,file = paste0(temp.dir,"p_value_file"),col.names = T,row.names = F,quote=F)
if(i1==1){
  res <- system(paste0("/data/zhangh24/software/plink2 --q-score-range ",temp.dir,"q_range_file ",temp.dir,"p_value_file header --threads 2 --score ",temp.dir,"prs_file header no-sum no-mean-imputation --bfile ",temp.dir,"all_chr --out ",temp.dir,"prs_",trait[l]))
  res = system(paste0("mv ",temp.dir,"*.profile ",out.dir))
}else{
  res <- system(paste0("/data/zhangh24/software/plink2 --q-score-range ",temp.dir,"q_range_file ",temp.dir,"p_value_file header --threads 2 --score ",temp.dir,"prs_file header no-sum no-mean-imputation --bfile ",temp.dir,"all_chr --out ",temp.dir,"prs_",trait[l],"_mega"))
  res = system(paste0("mv ",temp.dir,"*.profile ",out.dir))
}

