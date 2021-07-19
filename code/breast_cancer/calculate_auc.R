startend <- function(num,size,ind){
  split.all <- split(1:num,cut(1:num,size))
  temp <- split.all[[ind]]
  start <- temp[1]
  end <- temp[length(temp)]
  return(c(start,end))
}

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
#i1 for all SNPs
#i2 for mega SNPs
auc.result.vec = rep(0,2)

for(i1 in 1:2){
  
  
    pheno <- fread(paste0(data.dir,"GBHS_pheno.csv"))
    colnames(pheno)[5] = "ER_status"
    
    
    pheno = pheno %>% 
      mutate(overall_status=ifelse(Status=="Control",0,1)) %>% 
      filter(Age>=18&Age<99) %>% 
      mutate(Age = as.numeric(Age))
    set.seed(93)
    pheno = pheno[sample(c(1:nrow(pheno))),]
    
    
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
      #random order pheno.update so that cases and controls are split
      
      #cross-validation for the data
      pheno.update = left_join(pheno,prs.score,by=c("SB_ID"="ID"))
      for(i_rep in 1:n.rep){
        start.end <- startend(nrow(pheno.update),n.rep,i_rep)
        vad.id = c(start.end[1]:start.end[2])
        test.id = setdiff(c(1:nrow(pheno.update)),vad.id)
        test.data <- pheno.update[test.id,]
        vad.data <- pheno.update[vad.id,]
        roc_obj <- roc.binary(status="overall_status", variable=paste0("SCORE"),
                              confounders=~EV1+EV2+EV3+EV4+EV5+
                                EV6+EV7+EV8+EV9+EV10+Age,
                              data=test.data, precision=seq(0.1,0.9, by=0.1))
        auc.vec.test.rep[k,i_rep] = roc_obj$auc
        roc_obj <- roc.binary(status="overall_status", variable=paste0("SCORE"),
                              confounders=~EV1+EV2+EV3+EV4+EV5+
                                EV6+EV7+EV8+EV9+EV10+Age,
                              data=vad.data, precision=seq(0.1,0.9, by=0.1))
        
        auc.vec.vad.rep[k,i_rep] = roc_obj$auc
      }
      
    }
    
    #find the optimal C+T in each round
    best.vad.rep = rep(0,n.rep)
    for(i_rep in 1:n.rep){
      idx <- which.max(auc.vec.test.rep[,i_rep])
      best.vad.rep[i_rep] = auc.vec.vad.rep[idx,i_rep]
    }
    
    auc.result.vec[i1] = mean(best.vad.rep)  
  }
  
  
  


auc.result = auc.result.vec[2] 
auc.test.result = rowMeans(auc.vec.test.rep)
auc.vad.result = rowMeans(auc.vec.vad.rep)
auc.result = list(auc.result,auc.test.result,auc.vad.result)  
save(auc.result,file = paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/result/ct_auc_result_",l,"_mega"))


