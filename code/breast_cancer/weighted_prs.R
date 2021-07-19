startend <- function(num,size,ind){
  split.all <- split(1:num,cut(1:num,size))
  temp <- split.all[[ind]]
  start <- temp[1]
  end <- temp[length(temp)]
  return(c(start,end))
}
l = 1
setwd("/data/zhangh24/multi_ethnic/data/")
library(data.table)
library(tidyverse)
library(RISCA)
library(boot)
library(ROCnReg)


trait = c("overall","erpos","erneg","overall_ukb_meta")
#load best eur prs
prs.score = fread(paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/prs/prs_best_eur_",trait[l],".sscore"))
k =1 
prs.score = prs.score %>% 
  separate(IID,into=c("ID","ohter_id","nci_id"),remove=F)

prs.score.eur = prs.score %>% select(ID,paste0("SCORE",k,"_AVG")) %>% 
  rename(prs_eur = paste0("SCORE",k,"_AVG"))

#load best signle ethnic PRS
load(paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/result/ct_auc_result_",l,"_mega"))
idx <- which.max(auc.result[[2]])
out.dir.prs = "/data/zhangh24/multi_ethnic/result/breast_cancer/prs/"
prs.score = fread(paste0(out.dir.prs,"prs_",trait[l],"_mega.p_value_",idx,".profile"),header=T)
prs.score.tar = prs.score %>% 
  separate(IID,into=c("ID","ohter_id","nci_id"),remove=F) %>% 
  select(ID,SCORE) %>% 
  rename(prs_tar = SCORE)

data.dir = "/data/zhangh24/multi_ethnic/data/"
pheno <- fread(paste0(data.dir,"GBHS_pheno.csv"))
colnames(pheno)[5] = "ER_status"

pheno = pheno %>% 
  mutate(overall_status=ifelse(Status=="Control",0,1)) %>% 
  filter(Age>=18&Age<99) %>% 
  mutate(Age = as.numeric(Age))
set.seed(83)
pheno = pheno[sample(c(1:nrow(pheno))),]
n.rep = 2
auc.weighted.rep = rep(0,n.rep)
pheno.update = left_join(pheno,prs.score.eur,by=c("SB_ID"="ID")) %>% 
  left_join(prs.score.tar,by=c("SB_ID"="ID"))

for(i_rep in 1:n.rep){
  start.end <- startend(nrow(pheno.update),n.rep,i_rep)
  vad.id = c(start.end[1]:start.end[2])
  test.id = setdiff(c(1:nrow(pheno.update)),vad.id)
  test.data <- pheno.update[test.id,]
  vad.data <- pheno.update[vad.id,]
  
  model = glm(overall_status~prs_eur+prs_tar+EV1+EV2+EV3+EV4+EV5+
                EV6+EV7+EV8+EV9+EV10+Age,data=test.data)
  coef = coefficients(model)[2:3]
  prs_weight = cbind(vad.data$prs_eur,vad.data$prs_tar)%*%coef
  vad.data$prs_weight = prs_weight
  roc_obj <- roc.binary(status="overall_status", variable="prs_weight",
                        confounders=~EV1+EV2+EV3+EV4+EV5+
                          EV6+EV7+EV8+EV9+EV10+Age,
                        data=vad.data, precision=seq(0.1,0.9, by=0.1))
  
  auc.weighted.rep[i_rep] = roc_obj$auc
}

auc.weighted = mean(auc.weighted.rep)

model = glm(overall_status~prs_eur+prs_tar+EV1+EV2+EV3+EV4+EV5+
              EV6+EV7+EV8+EV9+EV10+Age,data=pheno.update)
coef = coefficients(model)[2:3]
prs_weight = cbind(pheno.update$prs_eur,pheno.update$prs_tar)%*%coef
pheno.update$prs_weight = prs_weight

roc_obj <- roc.binary(status="overall_status", variable="prs_weight",
                      confounders=~EV1+EV2+EV3+EV4+EV5+
                        EV6+EV7+EV8+EV9+EV10+Age,
                      data=pheno.update, precision=seq(0.1,0.9, by=0.1))
roc_obj$auc

