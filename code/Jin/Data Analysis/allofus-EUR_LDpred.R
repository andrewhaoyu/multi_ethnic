# ------------- Calculate cross-population PRS
rm(list=ls())
setwd('~/prs/')
library(data.table)
library(Rcpp)
library(readr)
library(MASS) # for mvrnorm
library(reshape) # for melt
library(parallel)
library(RcppTN)
library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(inline)
library(mvnfast)
library(genio) # for read_plink
library(dplyr)
library(stringr)
library(gdata)
library(R.utils) # for gzip
library(pROC)
library(bigsnpr)
library(readr)
maindir = '/dcs04/nilanjan/data/jjin/'
races = c('EUR','AFR','AMR')
traits = c('height','bmi')
temp <- commandArgs(TRUE)
trait = traits[as.numeric(temp[1])]
ldr = 3/1000
ncores = 1
set.seed(2020)
h2_seq <- c(0.7, 1, 1.4)
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
sets <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))

#r2 = pval = matrix(0,length(races), length(races))
for (trait in traits){
  for (chr in 1:22){
    for (j in 1:length(races)){
      race0 = races[j]
      scorefile = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race0,'/ldpred2/ldpred2-chr',chr,'.txt')
      if (file.exists(scorefile)){
        beta = bigreadr::fread2(scorefile)
        colnames(beta) = c(paste0('e',1:nrow(sets)),'a0','a1','rsid')
        rownames(beta) = beta$rsid
        beta0 = beta
        
        for (k in c(1:length(races))){
          race = races[k]
          tf = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race,'/ldpred2/crosspop/')
          if (!dir.exists(tf)) dir.create(tf)
          beta = beta0
          temfile = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race,'/ldpred2/crosspop/',race0,'-',race,'-chr',chr)
          # ------------- create score file
          bim = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jjin/UKBval/',race,'/genotype/chr',chr,'.bim'))
          bim = bim[!duplicated(bim[,2]),]
          # bfile
          bfile = paste0('/dcs04/nilanjan/data/jjin/UKBval/',race,'/genotype/chr',chr)
          a=bigreadr::fread2(paste0(bfile,'.bim'))
          dupid = a[duplicated(a[,2]),2]
          bim = bim[!(bim[,2]%in%dupid),]
          rownames(bim) = bim[,2]; colnames(bim) = c('chr','rsid','NA','pos','a0','a1')
          rsid = intersect(bim[,2], beta$rsid)
          beta = beta[rsid,]; bim = bim[rsid,]
          flipped = which((beta$a1 == bim$a0)&(beta$a0 == bim$a1));
          beta[flipped,c(1:nrow(sets))] = - beta[flipped,c(1:nrow(sets))];
          te = beta$a0; beta$a0[flipped] = beta$a1[flipped]; beta$a1[flipped] = te[flipped]; rm(te)
          prs.file = beta
          colnames(prs.file) = c(paste0('BETA',1:nrow(sets)),'REF', 'ALT', 'SNP')
          prs.file = prs.file[,c('SNP', 'ALT', paste0('BETA',1:nrow(sets)))]
          prs.file[is.na(prs.file)] = 0
          temdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race,'/ldpred2/crosspop/')
          if (!dir.exists(temdir)){dir.create(temdir)}
          prsscore = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race,'/ldpred2/crosspop/score-',race0,'-',race,'-chr',chr,'.txt')
          write.table(prs.file,file = prsscore,
                      col.names = T,row.names = F,quote=F)
          
          prscode = paste(paste0('/dcl01/chatterj/data/jin/software/plink2'),
                          paste0('--score ', prsscore),
                          'cols=+scoresums,-scoreavgs',
                          paste0('--score-col-nums 3-',nrow(sets)+2),
                          paste0('--bfile ',bfile),
                          " --threads 1",
                          " --allow-extra-chr",
                          paste0(' --out ', temfile))
          system(prscode)
          #}
        }
        rm(beta)
      }
    }
  }
}


# ------- Sum up across chromosomes
for (trait in traits){
  for (j in 1:length(races)){
    race0 = races[j]
    for (k in c(1:length(races))){
      race = races[k]
      #---------------------------------------#---------------------------------------
      chr = 22
      temfile = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race,'/ldpred2/crosspop/',race0,'-',race,'-chr',chr,'.sscore')
      if(file.exists(temfile)){
        y = bigreadr::fread2(temfile) 
        y[,paste0('SCORE',1:nrow(sets),'_SUM')] = y[,paste0('SCORE',1:nrow(sets),'_SUM')] # * n.snp[chr] # * dftem$ALLELE_CT # 
        print(paste0('Chr ', chr,' Completed'))
      }
      if(!file.exists(temfile)) print(chr)
      for(chr in c(21:1)){
        temfile = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race,'/ldpred2/crosspop/',race0,'-',race,'-chr',chr,'.sscore')
        if(file.exists(temfile)){
          dftem = bigreadr::fread2(temfile) 
          y[,paste0('SCORE',1:nrow(sets),'_SUM')] = y[,paste0('SCORE',1:nrow(sets),'_SUM')] + dftem[,paste0('SCORE',1:nrow(sets),'_SUM')] # * n.snp[chr] # * dftem$ALLELE_CT # 
          print(paste0('Chr ', chr,' Completed'))
        }
        if(!file.exists(temfile)) print(chr)
      }
      df.prs = y[,c('IID',paste0('SCORE',1:nrow(sets),'_SUM'))]
      colnames(df.prs) = c('id', paste0('prs',1:nrow(sets)))
      write.table(df.prs, paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race,'/ldpred2/crosspop/',race0,'-',race,'.txt'), 
                  row.names = F,col.names = T, quote = FALSE, sep = "\t" )
      rm(dftem)
      print(paste0('Complete ',race0,' on ',race))
    }
  }
}



# ------------------ Validation
traits = c('height','bmi')
races = c('EUR','AFR','AMR'); K = length(races)

Out = matrix(NA,length(traits),3*length(races))
colnames(Out) = sapply(races,function(x){paste(c('R2','SD','P-value'),x)})
rownames(Out) = traits
n.split = 1

h2_seq <- c(0.7, 1, 1.4) 
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
sets <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))

pars = matrix(NA,length(traits)*length(races),ncol(sets))
colnames(pars) = colnames(sets)
te = expand.grid(trait=traits,race=races)
rownames(pars) = sapply(1:nrow(te), function(x){paste(te[x,1], te[x,2])})
pars = as.data.frame(pars)



out.EUR = list()
load(paste0("/dcs04/nilanjan/data/jjin/prs/realdata/allofus/r2/r2-ldpred2.RData")) # Indx
R2 = matrix(NA, length(traits), K); colnames(R2) = races; rownames(R2) = traits
for (trait in traits){
  out = matrix(NA,n.split,3*length(races))
  colnames(out) = sapply(races,function(x){paste(c('R2 Adjusted','Regression Coef','P-value'),x)})
  for (race in races){
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
    temfile = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race,'/ldpred2/crosspop/EUR-',race,'.txt')
    preds = bigreadr::fread2(temfile)
    preds$id = as.character(preds$id)
    rownames(preds) = preds$id
    
    validatetable = validatetable[validatetable$id %in% rownames(preds),]
    
    valdat = validatetable; rm(validatetable)
    rownames(valdat) = as.character(valdat$id)
    #---------------------------------------#---------------------------------------
    
    # Validation
    validatetable = valdat[ids2,]
    
    output2 = matrix(rep(0,3),1,3)
    colnames(output2) = c('R2 Adjusted','Regression Coef','P-value')
    tem = data.frame(id=rownames(preds),prs = preds[,paste0('prs',Indx[trait,'EUR'])])
    tem$id = as.character(tem$id)
    rownames(tem) = tem$id
    tem = tem[validatetable$id,]
    
    prstable = merge(tem,validatetable,by='id')
    prstable = prstable[complete.cases(prstable),]
    if ((nrow(prstable)>0)&(sum(prstable$prs) != 0)){
      formula.prs=formula(paste0(paste('y', paste(c('age','gender',paste0('pc',1:10)),collapse="+"), sep='~')))
      fit = lm(formula.prs, data=prstable)
      prstable$res = fit$residuals
      fit = lm(res~prs, data=prstable)
      R2[trait,race] = summary(fit)$r.squared 
    }
  }
}

####### -------------------------------  tuned parameter:
save(R2, file=paste0("/dcs04/nilanjan/data/jjin/prs/realdata/allofus/r2/r2-eurldpred2.RData"))




