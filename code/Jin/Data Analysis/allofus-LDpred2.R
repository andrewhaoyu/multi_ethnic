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
library(bigparallelr)
maindir = '/dcs04/nilanjan/data/jjin/'
races = c('EUR','AFR','AMR')
traits = c('height','bmi')
temp <- commandArgs(TRUE)
race = races[as.numeric(temp[1])]
trait = traits[as.numeric(temp[2])]
chr =  as.numeric(temp[3])
ldr = 3/1000
ncores = 1


temdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait)
if (!dir.exists(temdir)){dir.create(temdir)}
temdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race)
if (!dir.exists(temdir)){dir.create(temdir)}

tedir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race,'/intermediate/')
if (!dir.exists(tedir)){dir.create(tedir)}

temdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race,'/ldpred2/')
if (!dir.exists(temdir)){dir.create(temdir)}
temfile = paste0(temdir,"ldpred2-chr",chr,".txt")
outfile = temfile
sumraw = bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/allofus_cleaned/",race,"/",trait,".txt"))
sumraw = sumraw[,c('rsID','CHR','POS_b38','N','BETA','SE','P','A1','A2')]
colnames(sumraw) = c('SNP_ID','CHR','POS','N','BETA','SE','PVAL','REF','ALT') # REF: allele corresponding to BETA

valdat = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jjin/UKBval/',race,'/genotype/chr',chr,'.bim'))
sumraw = sumraw[sumraw$SNP_ID %in% valdat[,2],]

#### read in reference data
refdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/refgeno/',race)
if (!dir.exists(refdir)){dir.create(refdir)}
temfile = paste0(refdir,'/chr',chr,'.bk')
system(paste0('rm -rf ',temfile))
temfile = paste0(refdir,'/chr',chr,'.rds')
system(paste0('rm -rf ',temfile))
snp_readBed(paste0(refdir,'/chr',chr,'.bed'))
obj.bigSNP <- snp_attach(paste0(refdir,'/chr',chr,'.rds'))
map <- obj.bigSNP$map[-c(3)]
names(map) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")

G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
POS2 <- snp_asGeneticPos(CHR, POS, dir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race,'/intermediate/'), ncores = ncores)
NCORES <-  nb_cores()


sumstats = sumraw[sumraw$CHR == chr,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]
set.seed(2020)
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
sumstats = sumstats[sumstats$rsid %in% map$rsid,]
info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
rownames(info_snp) = info_snp$rsid
ind.chr <- which(info_snp$chr == chr)
df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]

corr0 <- snp_cor(G, ind.col = ind.chr2, ncores = 3,
                 infos.pos = POS2[ind.chr2], size =  ldr)
corr <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))

# Automatic model
ldsc <- snp_ldsc2(corr0, df_beta)
h2_est <- ldsc[["h2"]]
print(paste0('Complete data preparation'))

h2_seq <- c(0.7, 1, 1.4)
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
params <- expand.grid(p = p_seq, h2 = signif(abs(h2_est) * h2_seq, 3), sparse = c(FALSE))

beta_grid <- snp_ldpred2_grid(corr, df_beta, params, burn_in = 50, num_iter = 200, ncores = NCORES)
rownames(beta_grid) = info_snp$rsid
beta_grid = cbind(beta_grid, info_snp[,c('a0','a1','rsid')])
colnames(beta_grid) = c(paste0('e',1:nrow(params)), 'a0','a1','rsid')

beta_grid[is.na(beta_grid)] = 0
beta_grid = as.data.frame(beta_grid)
write_delim(beta_grid,file = outfile, delim='\t')
rm(corr0, corr)
print(paste0('Complete'))



tedir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race,'/ldpred2/effect/')
if (!dir.exists(tedir)) dir.create(tedir)
tedir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race,'/ldpred2/score/')
if (!dir.exists(tedir)) dir.create(tedir)
outdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race,'/ldpred2/')

outfile = paste0(outdir, 'score/ldpred2-chr', chr,'.sscore')

# -------- PRS:
prs.file = data.frame(SNP = beta_grid$rsid, ALT = beta_grid$a0, BETA = beta_grid[,1:(ncol(beta_grid)-3)])
tem = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jjin/UKBval/',race,'/genotype/chr',chr,'.bim'))
dupid = tem[duplicated(tem[,2]),2]
prs.file = prs.file[!(prs.file$SNP %in% dupid),]
write.table(prs.file,file = paste0(outdir,'effect/ldpred2-chr',chr,'.txt'),
            col.names = T,row.names = F,quote=F)
prscode = paste(paste0('/dcl01/chatterj/data/jin/software/plink2'),
                paste0('--score ',  outdir,'effect/ldpred2-chr',chr,'.txt'),
                'cols=+scoresums,-scoreavgs',
                paste0('--score-col-nums 3-',ncol(prs.file)),
                paste0(' --bfile /dcs04/nilanjan/data/jjin/UKBval/',race,'/genotype/chr',chr),
                " --threads 1",
                paste0(' --out ', outdir, 'score/ldpred2-chr', chr))
print(paste0('Completed ',trait,race,chr))
system(prscode)

print(paste0('Complete Calculating Scores'))


# --------------------- Validation ---------------------
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
#library(rms)
library(DescTools)
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



Indx = data.frame(matrix(NA,length(traits),length(races)))
rownames(Indx) = traits; colnames(Indx) = races
R2 = matrix(NA, length(traits), K); colnames(R2) = races; rownames(R2) = traits
for (trait in traits){
  out = matrix(NA,n.split,3*length(races))
  colnames(out) = sapply(races,function(x){paste(c('R2 Adjusted','Regression Coef','P-value'),x)})
  for (race in races){
    writescore_path = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/allofus/',trait,'/',race,'/ldpred2/')
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
    temfile = paste0(writescore_path,'score/ldpred2.txt')
    preds = bigreadr::fread2(temfile)
    preds = preds[,c('IID',paste0('SCORE',1:nrow(sets),'_SUM'))]
    colnames(preds) = c('id',paste0('prs',1:nrow(sets)))
    preds$id = as.character(preds$id)
    rownames(preds) = preds$id
    
    validatetable = validatetable[validatetable$id %in% rownames(preds),]
    
    valdat = validatetable; rm(validatetable)
    rownames(valdat) = as.character(valdat$id)
    #---------------------------------------#---------------------------------------
    validatetable = valdat[ids1,]
    
    output = matrix(0,nrow(sets),3)
    colnames(output) = c('R2 Adjusted','Regression Coef','P-value')
    rownames(output) = sapply(1:nrow(sets),function(x){paste0('p=',sets[x,1],' h2=', sets[x,2], 'sp=',sets[x,3])})
    
    for(i in 1:nrow(sets)){
      tem = data.frame(id=rownames(preds), prs = preds[,paste0('prs',i)])
      prstable = merge(tem,validatetable,by='id')
      prstable = prstable[complete.cases(prstable),]
      if ((nrow(prstable)>0)&(sum(prstable$prs) != 0)){
        # get residual:
        formula.prs=formula(paste0(paste('y', paste(c('age','gender',paste0('pc',1:10)),collapse="+"), sep='~')))
        fit = lm(formula.prs, data=prstable)
        prstable$res = fit$residuals
        fit = lm(res~prs, data=prstable)
        output[i,'Regression Coef'] = coefficients(fit)['prs']
        output[i,'R2 Adjusted'] = summary(fit)$r.squared #(coefficients(fit)['prs'])^2/var(na.omit(prstable[,'res']))
        output[i,'P-value'] = summary(fit)$coefficients['prs','Pr(>|t|)']
        #print(i)
      }
    }
    Indx[trait,race] = which.max(output[,'R2 Adjusted'])
    tem = data.frame(id=rownames(preds), prs = preds[,paste0('prs',Indx[trait,race])])
    prstable = merge(tem,validatetable,by='id')
    prstable = prstable[complete.cases(prstable),]
    formula.prs=formula(paste0(paste('y', paste(c('age','gender',paste0('pc',1:10)),collapse="+"), sep='~')))
    fit = lm(formula.prs, data=prstable)
    prstable$res = fit$residuals
    fit.train = lm(res~prs, data=prstable)
    print(output[Indx[trait,race],])
    
    # Validation
    validatetable = valdat[ids2,]
    
    output2 = matrix(rep(0,3),1,3)
    colnames(output2) = c('R2 Adjusted','Regression Coef','P-value')
    tem = data.frame(id=rownames(preds),prs = preds[,paste0('prs',Indx[trait,race])])
    tem$id = as.character(tem$id)
    rownames(tem) = tem$id
    tem = tem[validatetable$id,]
    
    prstable = merge(tem,validatetable,by='id')
    prstable = prstable[complete.cases(prstable),]
    if ((nrow(prstable)>0)&(sum(prstable$prs) != 0)){
      formula.prs=formula(paste0(paste('y', paste(c('age','gender',paste0('pc',1:10)),collapse="+"), sep='~')))
      fit = lm(formula.prs, data=prstable)
      prstable$res = fit$residuals
      prstable$pred.prs = predict(fit.train, data.frame(prs = prstable$prs))
      fit = lm(res~pred.prs, data=prstable)
      output2[1,'Regression Coef'] = coefficients(fit)['pred.prs']
      output2[1,'R2 Adjusted'] = summary(fit)$r.squared # (coefficients(fit)['prs'])^2/var(na.omit(prstable[,'res']))
      output2[1,'P-value'] = summary(fit)$coefficients['pred.prs','Pr(>|t|)']
    }
    out[1,paste(c('R2 Adjusted','Regression Coef','P-value'),race)] = output2[1,c('R2 Adjusted','Regression Coef','P-value')]
    
    Out[trait,paste('R2',race)] = mean(out[,paste(c('R2 Adjusted'),race)])
    Out[trait,paste('SD',race)] = sqrt(var(out[,paste(c('R2 Adjusted'),race)]))
    Out[trait,paste('P-value',race)] = mean(out[,paste(c('P-value'),race)])
    print(Out[trait,])
    print(trait)
    pars[paste(trait,race),] = sets[Indx[trait,race],]
    R2[trait,race] = output2[1,'R2 Adjusted']
    save(pars,file=paste0("/dcs04/nilanjan/data/jjin/prs/realdata/allofus/r2/tuned-pars-ldpred2-",race,"-",trait,".RData"))
  }
}

####### -------------------------------  tuned parameter:
save(R2, Out, Indx, file=paste0("/dcs04/nilanjan/data/jjin/prs/realdata/allofus/r2/r2-ldpred2.RData"))




