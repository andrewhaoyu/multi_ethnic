startend <- function(num,size,ind){
  split.all <- split(1:num,cut(1:num,size))
  temp <- split.all[[ind]]
  start <- temp[1]
  end <- temp[length(temp)]
  return(c(start,end))
}
library(tidyverse)
library(data.table)
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
setwd("/data/zhangh24/multi_ethnic/data/")
library(data.table)
library(tidyverse)
library(RISCA)
library(boot)
library(caret)
library(SuperLearner)
library(ranger)
library(glmnet)

trait = c("overall","erpos","erneg")
data.dir = "/data/zhangh24/multi_ethnic/data/"
l = 1
out.dir.prs = "/data/zhangh24/multi_ethnic/result/breast_cancer/prs/"
n.rep = 2
#i1 for all SNPs
#i2 for mega SNPs

r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)

total = length(r2_vec)*length(wc_base_vec)*length(pthres)^2
  auc.vec.test.rep = matrix(0,total,n.rep)
  auc.vec.vad.rep = matrix(0,total,n.rep)
  
  pheno <- fread(paste0(data.dir,"GBHS_pheno.csv"))
  colnames(pheno)[5] = "ER_status"
  set.seed(83)
  #random order pheno.update so that cases and controls are split
  pheno = pheno[sample(c(1:nrow(pheno))),]
  pheno = pheno %>% 
    mutate(overall_status=ifelse(Status=="Control",0,1)) %>% 
    filter(Age>=18&Age<99) %>% 
    mutate(Age = as.numeric(Age))
  
  prs.mat = matrix(0,nrow(pheno),total)
  p.eur.vec = rep(0,total)
  p.tar.vec = rep(0,total)
  r2.vec = rep(0,total)
  wc.vec = rep(0,total)
  
  #pheno.update = pheno.update[sample(c(1:nrow(pheno.update))),]
  
  temp =1
  for(r_ind in 1:length(r2_vec)){
    wc_vec = wc_base_vec/r2_vec[r_ind]
    for(w_ind in 1:length(wc_vec)){
      for(k1 in 1:length(pthres)){
        for(k2 in 1:length(pthres)){
          prs.score = fread(paste0(out.dir.prs,"prs_",trait[l],"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,".p_value_",k2,".profile"),header=T)  
          prs.score = prs.score %>% 
            separate(IID,into=c("ID","ohter_id","nci_id"),remove=F)
          pheno.update = left_join(pheno,prs.score,by=c("SB_ID"="ID"))
          prs.mat[,temp] = pheno.update$SCORE
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
            auc.vec.test.rep[temp,i_rep] = roc_obj$auc
            roc_obj <- roc.binary(status="overall_status", variable=paste0("SCORE"),
                                  confounders=~EV1+EV2+EV3+EV4+EV5+
                                    EV6+EV7+EV8+EV9+EV10+Age,
                                  data=vad.data, precision=seq(0.1,0.9, by=0.1))
            
            auc.vec.vad.rep[temp,i_rep] = roc_obj$auc
          }
          p.eur.vec[temp] = k1
          p.tar.vec[temp] = k2
          r2.vec[temp] = r_ind
          wc.vec[temp] = w_ind
          temp = temp+1
          
          }
      }
    }
  }
  
    
    
  #find the optimal TDLD in each round
  best.vad.rep = rep(0,n.rep)
  for(i_rep in 1:n.rep){
    idx <- which.max(auc.vec.test.rep[,i_rep])
    best.vad.rep[i_rep] = auc.vec.vad.rep[idx,i_rep]
  }
  
  TDLD.result = mean(best.vad.rep)
  
  auc.test.vec = rowMeans(auc.vec.test.rep)
  auc.vad.vec = rowMeans(auc.vec.vad.rep)
  
  result.data <- data.frame(auc.test.vec,
                            auc.vad.vec,
                            p.eur.vec,
                            p.tar.vec,
                            r2.vec,
                            wc.vec)
  
  # SL.libray <- c(
  #   "SL.glmnet",
  #   #"SL.ridge",
  #   "SL.nnet"
  # )
  # prs.mat = as.data.frame(prs.mat)
  # prs.mat = prs.mat %>%
  #   select(which(!colSums(prs.mat)%in% 0))
  # mtx = cor(prs.mat)
  # drop = findCorrelation(mtx,cutoff=0.98)
  # drop = names(prs.mat)[drop]
  # prs.mat.new = prs.mat %>%
  #   select(-all_of(drop))
  # #prs.mat <- data.frame(ID = pheno.update$SB_ID,prs.mat.new)
  # auc.sl.rep = rep(0,n.rep)
  # for(i_rep in 1:n.rep){
  #   start.end <- startend(nrow(pheno.update),n.rep,i_rep)
  #   vad.id = c(start.end[1]:start.end[2])
  #   test.id = setdiff(c(1:nrow(pheno.update)),vad.id)
  #   test.data <- pheno.update[test.id,]
  #   vad.data <- pheno.update[vad.id,]
  #   y.test = test.data$overall_status
  #   y.vad = vad.data$overall_status
  #   x.test = prs.mat.new[test.id,]
  #   x.vad = prs.mat.new[vad.id,]
  #   sl = SuperLearner(Y = y.test, X = x.test, family = binomial(),
  #                     # For a real analysis we would use V = 10.
  #                     # V = 3,
  #                     SL.library = SL.libray)
  # 
  #   y.pred <- predict(sl,x.vad,onlySL = TRUE)
  # 
  #   vad.data = data.frame(vad.data,sl.prt = y.pred[[1]])
  # 
  #   roc_obj <- roc.binary(status="overall_status", variable=paste0("sl.prt"),
  #                         confounders=~EV1+EV2+EV3+EV4+EV5+
  #                           EV6+EV7+EV8+EV9+EV10+Age,
  #                         data=vad.data, precision=seq(0.1,0.9, by=0.1))
  #   auc.sl.rep[i_rep] = roc_obj$auc
  # }
  # 
auc.tdld = list(TDLD.result,result.data)
  

save(auc.tdld,file = paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/result/auc_tdld.rdata"))
