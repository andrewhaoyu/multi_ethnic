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
jindir = "/dcl01/chatterj/data/jin/"
RACES = c('EUR', 'AFR', 'AMR', 'EAS', 'SAS')
temp <- commandArgs(TRUE)
rho = as.numeric(temp[1])
GA = as.numeric(temp[2])
rep = as.numeric(temp[3])
race0 = 'EUR'; size0 = 4
temdir = paste0('/dcl01/chatterj/data/jin/prs/simulation/results/eurldpred/')
if (!dir.exists(temdir)){dir.create(temdir)}

# -------------------- only calculate race1 score on race2
for (chr in 1:22){
  tem = bigreadr::fread2(paste0('/dcl01/chatterj/data/jin/prs/simulation/results/ldpred2/ldpred2effect-mega-',race0,'-rho=',rho,
                                '-size=',size0,'-chr',chr,'-rep',rep,'-GA',GA,'.txt'))
  h2_seq <- 1:3
  p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
  params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))
  
  load(paste0("/dcs04/nilanjan/data/jjin/prs/sim/ldpred2/R2.ldpred2.RData"))
  pars = R2.ldpred2[(R2.ldpred2$rho == rho)&(R2.ldpred2$size == size0)&
                           (R2.ldpred2$GA == GA)&(R2.ldpred2$race == 'EUR'),]
  # ------------------ calculate PRS using PLINK ------------------
  col = which((params$p == pars$p)&(params$h2 == pars$h2))
  prs.file = data.frame(SNP = tem$marker.ID, ALT = tem$a0, BETA = tem[,paste0('e',col)])
  write.table(prs.file,file = paste0(temdir,'eurscore-',race0,'-rho=',rho,
                                     '-size=',size0,'-chr',chr,'-rep',rep,'-GA',GA,'.txt'),
              col.names = T,row.names = F,quote=F)
  for (race in RACES[-1]){
    for (size in 1:3){
      temf = paste0(temdir, 'eurprs-',race,'-rho=',rho,'-size=',size,'-chr',chr,'-rep',rep,'-GA',GA,'.sscore')
      prscode = paste(paste0('/dcl01/chatterj/data/jin/software/plink2'),
                      paste0('--score ', temdir, 'eurscore-',race0,'-rho=',rho,'-size=',size0,'-chr',chr,'-rep',rep,'-GA',GA,'.txt'),
                      paste0(' 1 2 3 --bfile /dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/chr',chr),
                      paste0('--threads 1'),
                      paste0('--allow-extra-chr'),
                      paste0('--keep /dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/IDs.eurtest-size',size,'.txt'),
                      paste0('--out ', temdir, 'eurprs-',race,'-rho=',rho,'-size=',size,'-chr',chr,'-rep',rep,'-GA',GA))
      system(prscode)
    }
  }
}




# ------------------ Validation: run locally ------------------
jindir = "/dcl01/chatterj/data/jin/"
RACES = c('EUR', 'AFR', 'AMR', 'EAS', 'SAS')[-1]
maindirec = "/dcl01/chatterj/data/yzhang/"
jindir = "/dcl01/chatterj/data/jin/"
settings = expand.grid(rho = 1:3, size = 1:4, GA = 1:5, race = RACES)
R2.ldpred2.EUR = cbind(settings,0,0,0,0,0)
colnames(R2.ldpred2.EUR)[5:9] = c('R2','P-value','p','p.common','h2')

size0 = 4
race0 = 'EUR'

# tuned LDpred parameters
load(paste0("/dcs04/nilanjan/data/jjin/prs/sim/ldpred2/R2.ldpred2.RData"))

for (set in 1:nrow(settings)){
  race = as.character(settings[set,'race'])
  rho = settings[set,'rho']
  size = settings[set,'size']
  GA = settings[set,'GA']
  # old: wrong (accidentally multiplied by n.SNP)
  pars = R2.ldpred2.mega[(R2.ldpred2.mega$rho == rho)&(R2.ldpred2.mega$size == size0)&
                           (R2.ldpred2.mega$GA == GA)&(R2.ldpred2.mega$race == 'EUR'),]
  for (rep in 1:10){
    eurfile = paste0('/dcl01/chatterj/data/jin/prs/simulation/results/eurldpred/eurprs-',race,'-rho=',rho,'-size=',size,'-rep',rep,'-GA',GA,'.txt')
    if (file.exists(eurfile)){
      eurprs = bigreadr::fread2(eurfile)
      if (sum(is.na(eurprs) == 0)){
        colnames(eurprs)[which(colnames(eurprs) == 'prs')] = 'eurprs'
        rownames(eurprs) = eurprs$id
        validatetable = read.table(paste0('/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/',race,'/pheno_summary_out_GA/phenotypes_rho',rho,'_',GA,'.phen'),header=F,
                                   colClasses = c(rep('character',2),rep('numeric',rep),rep('NULL',100-rep)))
        validatetable = as.data.frame(validatetable[,c(1,rep+2)])
        colnames(validatetable) = c('id','y')
        i=1
        dat = merge(validatetable, eurprs, by = 'id'); 
        dat$eurprs = scale(dat$eurprs,center=T,scale=T)
        
        # validation
        formula.prs=formula(paste0(paste(paste0('y'), paste(c('eurprs'),collapse="+"), sep='~')))
        fit = lm(formula.prs, data=dat)
        R2.ldpred2.EUR[set,c('R2','P-value','p','h2')] = R2.ldpred2.EUR[set,c('R2','P-value','p','h2')] + c((coefficients(fit)['eurprs'])^2/var(na.omit(dat[,paste0('y')])),
                                                                                                            summary(fit)$coefficients['eurprs','Pr(>|t|)'],
                                                                                                            pars[,c('p','h2')])
        print(set)
        print(R2.ldpred2.EUR[set,c('R2','P-value','p','h2')])
      }
    }
  }
}
R2.ldpred2.EUR = R2.ldpred2.EUR/10
save(R2.ldpred2.EUR,file=paste0("/dcl01/chatterj/data/jin/prs/simulation/results/eurldpred/R2.ldpred2.EUR.RData"))

