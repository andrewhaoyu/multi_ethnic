rm(list=ls())
library(readr)
library(data.table)
library(mvtnorm)
library(devtools)
library(lavaan)
library(gdata)
library(xtable)
library(MASS) # for the ginv
library(data.table)
library(corpcor) #for pseudoinverse
library(parallel)
library(MendelianRandomization) # for mr_ivw
library(dplyr)
library(R.utils) # for gzip
library(stringr) # for str_detect
library(genio) # a package to facilitate reading and writing genetics data. The focus of this vignette is processing plink BED/BIM/FAM files.
library(data.table)
library(pROC)
library(bigsnpr)
library(DescTools)
traits = c('HDL','LDL','TC','logTG')
races = c('EUR','AFR','AMR','EAS','SAS'); K = length(races)

out = matrix(0,length(races),length(races)+2)
colnames(out) = c('R2','P',paste0('Coef wprs',1:K))
rownames(out) = races

h2_seq <- c(0.7, 1, 1.4) 
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
sets <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))


# ------- Validation
out.wprs = matrix(NA,length(traits),K+1)
colnames(out.wprs) = c('R2 Adjusted',races); rownames(out.wprs) = traits
Indx = matrix(NA, length(traits), K^2)
colnames(Indx) = sapply(races, function(x){paste0(x,'-',races)})
rownames(Indx) = traits
R2 = matrix(NA, length(traits), K); colnames(R2) = races; rownames(R2) = traits
for (trait in traits){
  for (race in races[c(1:K)]){
    validatetable1 = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/',trait,'/tuning+validation/',race,'_tuning.txt'),header=T)
    ids1 = as.character(validatetable1$f.eid)
    validatetable2 = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/',trait,'/tuning+validation/',race,'_validation.txt'),header=T)
    ids2 = as.character(validatetable2$f.eid)
    validatetable = rbind(validatetable1[,1:2], validatetable2[,1:2])
    colnames(validatetable) = c('id','y')
    covariates1 = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/',race,'_tuning.txt'))
    covariates2 = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/',race,'_validation.txt'))
    covariates = rbind(covariates1, covariates2)
    colnames(covariates) = c('id','gender','age',paste0('pc',1:10))
    validatetable = merge(validatetable, covariates, by = 'id')
    validatetable = validatetable[complete.cases(validatetable$y),]
    validatetable$id = as.character(validatetable$id)
    valdat = validatetable
    
    rownames(valdat) = as.character(valdat$id)
    
    validatetable = valdat[valdat$id %in% ids1,]
    
    for (j in 1:K){ # training ancestry
      race0 = races[j]
      # tuning
      output = matrix(NA,nrow(sets),3)
      colnames(output) = c('R2','Coef','P')
      temprsfile = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait,'/',race,'/ldpred2/crosspop/',race0,'-',race,'.txt')
      tem = bigreadr::fread2(temprsfile)
      prstable = merge(tem,validatetable,by='id')
      
      prstable = prstable[complete.cases(prstable),]
      formula.prs=formula(paste0(paste('y', paste(c('age','gender',paste0('pc',1:10)),collapse="+"), sep='~')))
      fit = lm(formula.prs, data=prstable)
      prstable$res = fit$residuals
      for(i in 1:nrow(sets)){
        formula.prs = formula(paste0(paste('res', paste(c(paste0('prs',c(i))),collapse="+"), sep='~')))
        fit = lm(formula.prs, data=prstable)
        output[i,'R2'] = summary(fit)$r.squared #(coefficients(fit)[paste0('prs',i)])^2/var(na.omit(prstable[,'res']))
      }
      Indx[trait,paste0(race0,'-',race)] = which.max(output[,'R2'])
    }
    
    temtable = data.frame(id = validatetable$id)
    for (k in 1:K){
      race0 = races[k]
      temfile = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait,'/',race,'/ldpred2/crosspop/',race0,'-',race,'.txt')
      tem = bigreadr::fread2(temfile)
      tem = tem[,c('id', paste0('prs',Indx[trait,paste0(race0,'-',race)]))]
      rownames(tem) = tem$id
      colnames(tem) = c('id',paste0('prs',race0))
      tem = tem[validatetable$id,]
      temtable = merge(temtable, tem, by = 'id')
    }
    
    prstable = merge(temtable,validatetable,by='id')
    prstable = prstable[complete.cases(prstable),]
    if (nrow(prstable)>0){
      formula.prs=formula(paste0(paste('y', paste(c('age','gender',paste0('pc',1:10)),collapse="+"), sep='~')))
      fit.res = lm(formula.prs, data=prstable)
      prstable$res = fit.res$residuals
      formula.prs = formula(paste0(paste('res', paste(c(paste0('prs',races)),collapse="+"), sep='~')))
      fit.train = lm(formula.prs, data=prstable)
      for (race0 in races){
        out.wprs[trait,race0] = coefficients(fit.train)[paste0('prs',race0)]
      }
    }
    
    # validation
    validatetable = valdat[valdat$id %in% ids2,]
    temtable = data.frame(id = validatetable$id)
    for (k in 1:K){
      race0 = races[k]
      temfile = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait,'/',race,'/ldpred2/crosspop/',race0,'-',race,'.txt')
      tem = bigreadr::fread2(temfile)
      tem = tem[,c('id', paste0('prs',Indx[trait,paste0(race0,'-',race)]))]
      rownames(tem) = tem$id
      colnames(tem) = c('id',paste0('prs',race0))
      tem = tem[validatetable$id,]
      temtable = merge(temtable, tem, by = 'id')
    }
    prstable = merge(temtable,validatetable,by='id')
    prstable = prstable[complete.cases(prstable),]
    if (nrow(prstable)>0){
      formula.prs=formula(paste0(paste('y', paste(c('age','gender',paste0('pc',1:10)),collapse="+"), sep='~')))
      fit = lm(formula.prs, data=prstable)
      prstable$res = fit$residuals
      pred = as.matrix(prstable[,c(paste0('prs',races))]) %*% t(t(summary(fit.train)$coef[c(paste0('prs',races)),'Estimate']))
      prstable$prs = pred
      fit = lm(res~prs, data=prstable)
      R2[trait,race] = summary(fit)$r.squared
    }
    print(R2)
    print(paste0(trait,' completed.'))
  }
}

save(R2,file=paste0("/dcs04/nilanjan/data/jjin/prs/realdata/glgc/r2/r2-wprs.RData"))




