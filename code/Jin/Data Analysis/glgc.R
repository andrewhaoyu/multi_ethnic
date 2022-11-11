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
races = c('EUR','AFR','AMR','EAS','SAS')
traits = c('HDL','LDL','TC','logTG')
temp <- commandArgs(TRUE)
race = races[as.numeric(temp[1])]
trait = traits[as.numeric(temp[2])]
chr =  as.numeric(temp[3])
ldr = 3/1000
ncores = 1


temdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait)
if (!dir.exists(temdir)){dir.create(temdir)}
temdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait,'/',race)
if (!dir.exists(temdir)){dir.create(temdir)}
temdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait,'/',race,'/ldpred2/')
if (!dir.exists(temdir)){dir.create(temdir)}
temfile = paste0(temdir,"ldpred2-chr",chr,".txt")
outfile = temfile

topsnp = bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/hugeh2_locus/",trait,"_topSNP_corrected.txt"))
sumraw = bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",race,"/",trait,".txt"))
sumtop = bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",race,"/",trait,".txt"))
sumtop = sumtop[(sumtop$rsID %in% topsnp$topSNP),]
sumraw = sumraw[,c('rsID','CHR','POS_b37','N','BETA','SE','P','A1','A2')]
colnames(sumraw) = c('SNP_ID','CHR','POS','N','BETA','SE','PVAL','REF','ALT') # REF: allele corresponding to BETA
valdat = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jjin/UKBval/',race,'/genotype/chr',chr,'.bim'))
sumraw = sumraw[sumraw$SNP_ID %in% valdat[,2],]

#### read in reference data
temdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/refgeno/',race,'/intermediate/')
if (!dir.exists(temdir)){dir.create(temdir)}
temdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/refgeno/',race,'/intermediate/',trait)
if (!dir.exists(temdir)){dir.create(temdir)}

refdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/refgeno/',race)
if (!dir.exists(refdir)){dir.create(refdir)}
temfile = paste0(refdir,'/chr',chr,'.bk')
system(paste0('rm -rf ',temfile))
temfile = paste0(refdir,'/chr',chr,'.rds')
system(paste0('rm -rf ',temfile))
snp_readBed(paste0(refdir,'/chr',chr,'.bed'))
obj.bigSNP <- snp_attach(paste0(refdir,'/chr',chr,'.rds'))
map <- obj.bigSNP$map[-c(3)]
names(map) <- c("chr", "rsid", "pos", "a0", "a1")

G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
POS2 <- snp_asGeneticPos(CHR, POS, dir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/refgeno/',race,'/intermediate/',trait), ncores = ncores)
NCORES <-  nb_cores()


sumstats = sumraw[sumraw$CHR == chr,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]
set.seed(2020)
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
sumstats = sumstats[sumstats$rsid %in% map$rsid,]
info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
rownames(info_snp) = info_snp$rsid
## compute correlation
## indices in info_snp
ind.chr <- which(info_snp$chr == chr)
df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
## indices in G
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
colnames(beta_grid) = c(paste0('e',1:nrow(params)), 'a0','a1','rsid') # beta corresponds to a0!

beta_grid[is.na(beta_grid)] = 0 # LDpred2 sometimes gives NA output
write_delim(beta_grid,file = outfile, delim='\t')
rm(corr0, corr)
print(paste0('Complete'))





# ---------------------- Calculate cross-population PRS ----------------------
maindir = '/dcs04/nilanjan/data/jjin/'
traits = c('HDL','LDL','TC','logTG')
races = c('EUR','AFR','AMR','EAS','SAS')
temp <- commandArgs(TRUE)
trait = traits[as.numeric(temp[1])]
ldr = 3/1000
ncores = 1
set.seed(2020)
h2_seq <- c(0.7, 1, 1.4)
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
sets <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))

for (trait in traits){
  topsnp = bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/hugeh2_locus/",trait,"_topSNP_corrected.txt"))
  for (chr in 1:22){
    for (j in 1:length(races)){
      race0 = races[j]
      scorefile = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait,'/',race0,'/ldpred2/effect/ldpred2-chr', chr,'.txt')
      if (file.exists(scorefile)){
        beta = bigreadr::fread2(scorefile)
        colnames(beta) = c('rsid','a1',paste0('e',1:nrow(sets)))
        rownames(beta) = beta$rsid
        
        snpid = which(beta$rsid %in% topsnp$topSNP)
        beta0 = beta
        if (length(snpid) > 0){
          beta = beta0[-snpid,]
        }
        beta0 = beta
        
        for (k in c(1:length(races))){
          race = races[k]
          beta = beta0
          temfile = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait,'/',race,'/ldpred2/crosspop/',race0,'-',race,'-chr',chr)
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
          beta[flipped,c(3:ncol(beta))] = - beta[flipped,c(3:ncol(beta))];
          te = beta$a0; beta$a0[flipped] = beta$a1[flipped]; beta$a1[flipped] = te[flipped]; rm(te)
          prs.file = beta
          colnames(prs.file) = c('SNP', 'ALT', paste0('BETA',1:nrow(sets)))
          prs.file[is.na(prs.file)] = 0
          temdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait,'/',race,'/ldpred2/crosspop/')
          if (!dir.exists(temdir)){dir.create(temdir)}
          prsscore = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait,'/',race,'/ldpred2/crosspop/score-',race0,'-',race,'-chr',chr,'.txt')
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
        }
        rm(beta)
      }
    }
  }
}



for (trait in traits){
  topsnp = bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/hugeh2_locus/",trait,"_topSNP_corrected.txt"))
  for (j in 1:length(races)){
    race0 = races[j]
    for (k in c(1:length(races))){
      race = races[k]
      #---------------------------------------#---------------------------------------
      chr = 22
      temfile = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait,'/',race,'/ldpred2/crosspop/',race0,'-',race,'-chr',chr,'.sscore')
      if(file.exists(temfile)){
        y = bigreadr::fread2(temfile) 
        y[,paste0('SCORE',1:nrow(sets),'_SUM')] = y[,paste0('SCORE',1:nrow(sets),'_SUM')] # * n.snp[chr] # * dftem$ALLELE_CT # 
        print(paste0('Chr ', chr,' Completed'))
      }
      if(!file.exists(temfile)) print(chr)
      for(chr in c(21:1)){
        temfile = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait,'/',race,'/ldpred2/crosspop/',race0,'-',race,'-chr',chr,'.sscore')
        if(file.exists(temfile)){
          dftem = bigreadr::fread2(temfile) 
          y[,paste0('SCORE',1:nrow(sets),'_SUM')] = y[,paste0('SCORE',1:nrow(sets),'_SUM')] + dftem[,paste0('SCORE',1:nrow(sets),'_SUM')] # * n.snp[chr] # * dftem$ALLELE_CT # 
          print(paste0('Chr ', chr,' Completed'))
        }
        if(!file.exists(temfile)) print(chr)
      }
      df.prs = y[,c('IID',paste0('SCORE',1:nrow(sets),'_SUM'))]
      colnames(df.prs) = c('id', paste0('prs',1:nrow(sets)))
      write.table(df.prs, paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait,'/',race,'/ldpred2/crosspop/',race0,'-',race,'.txt'), 
                  row.names = F,col.names = T, quote = FALSE, sep = "\t" )
      rm(dftem)
      print(paste0('Complete ',race0,' on ',race))
    }
  }
}


for (trait in traits){
  topsnp = bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/hugeh2_locus/",trait,"_topSNP_corrected.txt"))
  outdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait,'/crosstopsnp/')
  for (j in 1:length(races)){
    race0 = races[j]
    scoredir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait,'/',race0,'/ldpred2/')
    sumtop = bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",race0,"/",trait,".txt"))
    sumtop = sumtop[(sumtop$rsID %in% topsnp$topSNP),]
    chrs = unique(sumtop$CHR)
    for (k in c(1:length(races))){
      race = races[k]
      
      chr = chrs[1]
      temfile = paste0(outdir, race0,'-',race,'-topsnp-chr',chr,'.sscore')
      if(file.exists(temfile)){
        y = bigreadr::fread2(temfile)
        y[is.na(y)] = 0
        print(paste0('Top SNP Chr ', chr,' Completed'))
      }
      if(!file.exists(temfile)) print(chr)
      for (chr in chrs){
        temfile = paste0(outdir, race0,'-',race,'-topsnp-chr',chr,'.sscore')
        if(file.exists(temfile)){
          pred_grid = bigreadr::fread2(temfile)
          pred_grid[is.na(pred_grid)] = 0
          y[,5:ncol(y)] = y[,5:ncol(y)] + pred_grid[,5:ncol(y)]
          print(paste0('Top SNP: Chr ', chr,' Completed'))
        }
        if(!file.exists(temfile)) print(chr)
      }
      write.table(y,paste0(outdir, race0,'-',race,'-topsnp.txt'), row.names = F,col.names = T, quote = FALSE, sep = "\t" )
      rm(y)
    }
  }
}





# ----------------- Validation, run locally -----------------
traits = c('HDL','LDL','TC','logTG')
races = c('EUR','AFR','AMR','EAS','SAS'); K = length(races)
maindirec = "/dcl01/chatterj/data/yzhang/"
jindir = "/dcl01/chatterj/data/jin/"

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
    writescore_path = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/',trait,'/',race,'/ldpred2/')
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
        prstable$prs = scale(prstable$prs,center=T,scale=T)
        # get residual:
        formula.prs=formula(paste0(paste('y', paste(c('age','gender',paste0('pc',1:10)),collapse="+"), sep='~')))
        fit = lm(formula.prs, data=prstable)
        prstable$res = fit$residuals
        fit = lm(res~prs, data=prstable)
        output[i,'Regression Coef'] = coefficients(fit)['prs']
        output[i,'R2 Adjusted'] = summary(fit)$r.squared #(coefficients(fit)['prs'])^2/var(na.omit(prstable[,'res']))
        output[i,'P-value'] = summary(fit)$coefficients['prs','Pr(>|t|)']
      }
    }
    Indx[trait,race] = which.max(output[,'R2 Adjusted']) #which.max(output[,'Regression Coef'])
    
    # Validation
    validatetable = valdat[ids2,]
    
    output2 = matrix(rep(0,3),1,3)
    colnames(output2) = c('R2 Adjusted','Regression Coef','P-value')
    tem = data.frame(id=rownames(preds),prs = preds[,paste0('prs',Indx[trait,race])])
    rownames(tem) = tem$id
    tem = tem[validatetable$id,]
    
    prstable = merge(tem,validatetable,by='id')
    prstable = prstable[complete.cases(prstable),]
    if ((nrow(prstable)>0)&(sum(prstable$prs) != 0)){
      prstable$prs = scale(prstable$prs,center=T,scale=T)
      formula.prs=formula(paste0(paste('y', paste(c('age','gender',paste0('pc',1:10)),collapse="+"), sep='~')))
      fit = lm(formula.prs, data=prstable)
      prstable$res = fit$residuals
      fit = lm(res~prs, data=prstable)
      R2[trait,race] = summary(fit)$r.squared # (coefficients(fit)['prs'])^2/var(na.omit(prstable[,'res']))
    }
    pars[paste(trait,race),] = sets[Indx[trait,race],]
    print(trait)
    print(R2)
    save(pars,file=paste0("/dcs04/nilanjan/data/jjin/prs/realdata/glgc/r2-notopsnp/tuned-pars-ldpred2-",race,"-",trait,".RData"))
  }
}
####### -------------------------------  tuned parameter:
save(R2, Out, Indx, file=paste0("/dcs04/nilanjan/data/jjin/prs/realdata/glgc/r2/summary-ldpred2.RData"))



