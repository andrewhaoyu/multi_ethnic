# --------------- Calculate cross-ancestry scores based on LDpred2 effect sizes from LDpred2.R:
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
maindir = '/dcs04/nilanjan/data/jjin/'
races = c('EUR', 'AFR', 'AMR', 'EAS', 'SAS')
temp <- commandArgs(TRUE)
chr =  as.numeric(temp[1])
rho = as.numeric(temp[2])
size = as.numeric(temp[3])
GA = as.numeric(temp[4])
rep = as.numeric(temp[5])
size0 = 4

ldr = 3/1000
ncores = 1
set.seed(2020)
h2_seq <- 1:3
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
sets <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))
load(paste0("/dcs04/nilanjan/data/jjin/prs/sim/ldpred2/R2.ldpred2.mega.rep",rep,".RData"))

temdir = paste0("/dcs04/nilanjan/data/jjin/prs/sim/crosspop/")
if (!dir.exists(temdir)){dir.create(temdir)}
temdir = paste0("/dcs04/nilanjan/data/jjin/prs/sim/crosspop/score/")
if (!dir.exists(temdir)){dir.create(temdir)}

for (j in 1:length(races)){
  race0 = races[j]
  if (race0 == 'EUR'){
    scorefile = paste0('/dcl01/chatterj/data/jin/prs/simulation/results/ldpred2/ldpred2effect-mega-',
                       race0,'-rho=',rho,'-size=',size0,'-chr',chr,'-rep',rep,'-GA',GA,'.txt')
  }
  if (race0 != 'EUR'){
    scorefile = paste0('/dcl01/chatterj/data/jin/prs/simulation/results/ldpred2/ldpred2effect-mega-',
                       race0,'-rho=',rho,'-size=',size,'-chr',chr,'-rep',rep,'-GA',GA,'.txt')
  }
  scores = bigreadr::fread2(scorefile)
  beta = scores
  colnames(beta)[1] = 'rsid'
  rownames(beta) = beta$rsid
  for (k in c(1:length(races))[-c(j)]){
    race = races[k]
    tfile = paste0(temdir,race0,'-',race,'-rho',rho,'-size',size,'-rep',rep,'-GA',GA,'-chr',chr,'.sscore')
    if (!file.exists(tfile)){
      #k=1
      testsize = ifelse(race == 'EUR', 4, size)
      # ------------- create score file
      bim = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jjin/prs/sim/',race,'/geno/mega/chr',chr,'.bim'))
      bim = bim[!duplicated(bim[,2]),] # same as below
      # bfile
      bfile = paste0('/dcs04/nilanjan/data/jjin/prs/sim/',race,'/geno/mega/chr',chr)
      a=bigreadr::fread2(paste0(bfile,'.bim'))
      dupid = a[duplicated(a[,2]),2]
      bim = bim[!(bim[,2]%in%dupid),]
      rownames(bim) = bim[,2]; colnames(bim) = c('chr','rsid','NA','pos','a0','a1')
      rsid = intersect(bim[,2], beta$rsid)
      beta = beta[rsid,]; bim = bim[rsid,]
      #matched = which((beta$a0 == bim$a0)&(beta$a1 == bim$a1)); 
      flipped = which((beta$a1 == bim$a0)&(beta$a0 == bim$a1));
      beta[flipped,4:ncol(beta)] = - beta[flipped,4:ncol(beta)];
      te = beta$a0; beta$a0[flipped] = beta$a1[flipped]; beta$a1[flipped] = te[flipped]; rm(te)
      prs.file = beta[,c(1,3,4:ncol(beta))]
      colnames(prs.file) = c('SNP','ALT',paste0('BETA',1:nrow(sets)))
      prsscore = paste0(temdir,race0,'-',race,'-rho',rho,'-size',size,'-rep',rep,'-GA',GA,'-chr',chr,'.txt')
      write.table(prs.file,file = prsscore,
                  col.names = T,row.names = F,quote=F)
      
      temfile = paste0(temdir,race0,'-',race,'-rho',rho,'-size',size,'-rep',rep,'-GA',GA,'-chr',chr)
      prscode = paste(paste0('/dcl01/chatterj/data/jin/software/plink2'),
                      paste0('--score ', prsscore),
                      'cols=+scoresums,-scoreavgs',
                      paste0('--score-col-nums 3-',nrow(sets)+2),
                      paste0('--bfile ',bfile),
                      " --threads 1",
                      " --allow-extra-chr",
                      paste0('--keep /dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/prs_size',testsize,'.id.txt'),
                      paste0(' --out ', temfile))
      system(prscode)
    }
  }
}


# ------------------------ Validation, run locallly ------------------------
jindir = "/dcl01/chatterj/data/jin/"
races = c('EUR','AFR','AMR','EAS','SAS')
maindirec = "/dcl01/chatterj/data/yzhang/"
jindir = "/dcl01/chatterj/data/jin/"
K = length(races)
R2.wldpred2 = matrix(rep(0,8),nrow=1)
colnames(R2.wldpred2) = c('R2','Regression Coef','P-value','p1','p2','p3','p4','r12')
h2_seq <- 1:3
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
sets <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))

n.training = c(15000, 45000, 80000, 100000)
n.total = 1.2e5
split = list()
for (ethnicity in races){
  split[[ethnicity]] = list()
  for (siz in 1:4){
    split[[ethnicity]][[siz]] = list()
    n.test = n.val = (n.total-n.training[siz])/2#n.training[siz]*0.1
    sample.test = (n.training[siz]+1):(n.training[siz]+n.val)
    sample.val = (n.training[siz]+n.val+1):(n.total)
    split[[ethnicity]][[siz]] = list(1:n.training[siz], sample.test, sample.val)
  }
}
rep = 1
settings = rbind(expand.grid(rho = 1:3, size = 1:4, GA = 1:5, race = races[-1]))
writescore_path = paste0('/dcs04/nilanjan/data/jjin/prs/sim/ldpred2/')

out = matrix(0,length(races),2*length(races)+2)
colnames(out) = c('R2','P',paste0('Coef wprs',1:K), paste0('P wprs',1:K))
rownames(out) = races

Indx = matrix(NA, nrow(settings), K)
colnames(Indx) = races
R2.wldpred2 = cbind(settings,0,0,0)
colnames(R2.wldpred2)[(ncol(R2.wldpred2)-2): (ncol(R2.wldpred2))] = c('R2','P-value','p')

for (set in 1:nrow(settings)){
  race = as.character(settings[set,'race'])
  rho = settings[set,'rho']
  size = settings[set,'size']
  size1 = 4; size2 = size
  sizes = c(size1,rep(size2,K-1))
  GA = settings[set,'GA']
  
  for (rep in 1:10){
    validatetable = read.table(paste0('/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/',race,'/pheno_summary_out_GA/phenotypes_rho',rho,'_',GA,'.phen'),header=F,
                               colClasses = c(rep('character',2),rep('numeric',rep),rep('NULL',100-rep)))
    validatetable = as.data.frame(validatetable[,c(1,rep+2)])
    colnames(validatetable) = c('id','y')
    valid = bigreadr::fread2(paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/test_size',sizes[which(races == race)],'.id.txt'))
    traindat = validatetable[valid[,1],] # should delete in the correct version
    
    for (j in 1:K){ # training ancestry
      race0 = races[j]
      # tuning
      output = matrix(NA,nrow(sets),3)
      colnames(output) = c('R2','Coef','P')
      if (race == race0){
        temprsfile = paste0(writescore_path,'mega-',race,'-rho=',rho,'-size=',size,'-rep',rep,'-GA',GA,'.txt')
        tem = bigreadr::fread2(temprsfile)
        colnames(tem) = paste0('prs',1:nrow(sets))
        tem$id = as.character(c(split[[race]][[size]][[2]],split[[race]][[size]][[2]]))
        tem[is.na(tem)] = 0
      } 
      if (race != race0){
        temprsfile = paste0('/dcs04/nilanjan/data/jjin/prs/sim/crosspop/score/prs-',race0,'-',race,'-rho',rho,'-size',size,'-rep',rep,'-GA',GA,'.txt')
        tem = bigreadr::fread2(temprsfile)
        colnames(tem)[6:ncol(tem)] = paste0('prs',1:nrow(sets))
      } 
      
      prstable = merge(tem,traindat,by='id')
      prstable = prstable[complete.cases(prstable),]
      for(i in 1:nrow(sets)){
        formula.prs = formula(paste0(paste('y', paste(c(paste0('prs',i)),collapse="+"), sep='~')))
        fit = lm(formula.prs, data=prstable)
        output[i,'R2'] = summary(fit)$r.squared #(coefficients(fit)[paste0('prs',i)])^2/var(na.omit(prstable[,'res']))
      }
      Indx[set,race0] = which.max(output[,'R2'])
    }
    
    temtable = data.frame(id = traindat$id)
    for (k in 1:K){
      race0 = races[k]
      if (race == race0){
        temprsfile = paste0(writescore_path,'mega-',race,'-rho=',rho,'-size=',size,'-rep',rep,'-GA',GA,'.txt')
        tem = bigreadr::fread2(temprsfile)
        colnames(tem) = paste0('prs',1:nrow(sets))
        tem$id = as.character(c(split[[race]][[size]][[2]],split[[race]][[size]][[3]]))
      } 
      if (race != race0){
        temprsfile = paste0('/dcs04/nilanjan/data/jjin/prs/sim/crosspop/score/prs-',race0,'-',race,'-rho',rho,'-size',size,'-rep',rep,'-GA',GA,'.txt')
        tem = bigreadr::fread2(temprsfile)
        colnames(tem)[6:ncol(tem)] = paste0('prs',1:nrow(sets))
      } 
      tem = tem[,c('id', paste0('prs',Indx[set,race0]))]
      tem$id = as.character(tem$id)
      rownames(tem) = tem$id
      colnames(tem) = c('id',paste0('prs',race0))
      tem = tem[traindat$id,]
      temtable = merge(temtable, tem, by = 'id')
    }
    
    output = matrix(NA,1,4+2)
    colnames(output) = c('R2', 'P', sapply(1:4, function(x){paste0(c('Coef prs'), x)}))
    
    prstable = merge(traindat,temtable,by='id');
    prstable = prstable[complete.cases(prstable),]
    formula.prs = formula(paste0(paste('y', paste(c(paste0('prs',races)),collapse="+"), sep='~')))
    fit = lm(formula.prs, data=prstable)
    output[1,paste0('R2')] = summary(fit)$r.squared #(coefficients(fit)['wprs'])^2/var(na.omit(prstable2[,'res']))
    output[1,paste0('P')] = pf(summary(fit)$fstatistic[1],df1 = summary(fit)$fstatistic[2], df2 = summary(fit)$fstatistic[3], lower.tail=F)
    fit.train = fit
    
    # -------- validation
    valid = bigreadr::fread2(paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/validation_size',sizes[which(races == race)],'.id.txt'))
    valdat = validatetable[valid[,1],]
    temtable = data.frame(id = valdat$id)
    for (k in 1:K){
      race0 = races[k]
      if (race == race0){
        temprsfile = paste0(writescore_path,'mega-',race,'-rho=',rho,'-size=',size,'-rep',rep,'-GA',GA,'.txt')
        tem = bigreadr::fread2(temprsfile)
        colnames(tem) = paste0('prs',1:nrow(sets))
        tem$id = as.character(c(split[[race]][[size]][[2]],split[[race]][[size]][[3]]))
      } 
      if (race != race0){
        temprsfile = paste0('/dcs04/nilanjan/data/jjin/prs/sim/crosspop/score/prs-',race0,'-',race,'-rho',rho,'-size',size,'-rep',rep,'-GA',GA,'.txt')
        tem = bigreadr::fread2(temprsfile)
        colnames(tem)[6:ncol(tem)] = paste0('prs',1:nrow(sets))
      } 
      tem = tem[,c('id', paste0('prs',Indx[set,race0]))]
      tem$id = as.character(tem$id)
      rownames(tem) = tem$id
      colnames(tem) = c('id',paste0('prs',race0))
      temtable = merge(temtable, tem, by = 'id')
    }
    
    prstable = merge(temtable,valdat,by='id');
    prstable$pred = predict(fit.train,prstable)
    fit = lm(y~pred,data=prstable)
    R2.wldpred2[set,c('R2','P-value','p')] = R2.wldpred2[set,c('R2','P-value','p')] + c(summary(fit)$r.squared, summary(fit)$coef['pred','Pr(>|t|)'], 0)
    print(set)
  }
}
R2.wldpred2 = R2.wldpred2/10
save(R2.wldpred2, Indx, file=paste0("/dcs04/nilanjan/data/jjin/prs/sim/r2new/R2.weighted_ldpred2.RData"))

