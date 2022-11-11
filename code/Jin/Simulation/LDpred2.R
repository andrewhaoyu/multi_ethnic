rm(list=ls())
setwd('/fastscratch/myscratch/jjin/')
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
temp <- commandArgs(TRUE)
RACES = c('EUR', 'AFR', 'AMR', 'EAS', 'SAS')
race = RACES[as.numeric(temp[1])]
chr =  as.numeric(temp[2])
rho = as.numeric(temp[3])
size = as.numeric(temp[4])
GA = as.numeric(temp[5])
rep = as.numeric(temp[6])

ldr = 3/1000
sum.raw = bigreadr::fread2(paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/sumdata/megasum-rho',rho,'-size',size,'-rep',rep,'-GA',GA,'-chr',chr,'.txt'))
sum.raw = as.data.frame(sum.raw)

# ------------------------ Run LDpred2
temfile = paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/ref_chr',chr,'.bk')
system(paste0('rm -rf ',temfile))
temfile = paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/ref_chr',chr,'.rds')
system(paste0('rm -rf ',temfile))
temfile = paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/ref_chr',chr,'.bk')
if (!file.exists(temfile)){
  snp_readBed(paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/ref_chr',chr,'.bed'))
}
obj.bigSNP <- snp_attach(paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/ref_chr',chr,'.rds'))

G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
NCORES <- nb_cores()
# Read external summary statistics
sumstats = sum.raw[sum.raw$CHR == chr,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]
#str(sumstats)

set.seed(2020)
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
map <- obj.bigSNP$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats, map, strand_flip = T)
rownames(info_snp) = info_snp$rsid
## compute correlation
POS2 <- snp_asGeneticPos(CHR, POS, dir = paste0('/dcs04/nilanjan/data/jjin/prs/sim/',race,'/sumdata/'), ncores = 2)
## indices in info_snp
ind.chr <- which(info_snp$chr == chr)
df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
## indices in G
ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
corr0 <- snp_cor(G, ind.col = ind.chr2, ncores = 4, #size = ldradius)
                 infos.pos = POS2[ind.chr2], size = ldr) # default
corr <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))

# Automatic model
ldsc <- snp_ldsc2(corr0, df_beta)
h2_est <- abs(ldsc[["h2"]])


# grid of models:
H2_seq = h2_seq <- signif(abs(h2_est) * c(0.7, 1, 1.4), 3)
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))

beta_grid <- snp_ldpred2_grid(corr, df_beta, params, burn_in = 50, num_iter = 200, ncores = NCORES)
beta_grid = as.data.frame(beta_grid)
rownames(beta_grid) = info_snp$rsid
beta_grid = cbind(info_snp$rsid, info_snp$a0, info_snp$a1, beta_grid)
colnames(beta_grid) = c(c('marker.ID', 'a0', 'a1'),paste0('e',1:nrow(params)))



beta_grid = beta_grid[,-2]
beta_grid[is.na(beta_grid)] = 0
# bigreadr::fwrite2(beta_grid, paste0('/dcs04/nilanjan/data/jjin/prs/test/ldpred2effect-mega-',race,'-rho=',rho,
#                                     '-size=',size,'-chr',chr,'-rep',rep,'-GA',GA,'.txt'), colnames = F)
write.table(beta_grid,file = paste0('/dcs04/nilanjan/data/jjin/prs/sim/results/ldpred2/ldpred2prs-mega-',race,'-rho=',rho,
                                    '-size=',size,'-chr',chr,'-rep',rep,'-GA',GA,'.txt'),
            col.names = T,row.names = F,quote=F)
prscode = paste(paste0('/dcl01/chatterj/data/jin/software/plink2'),
                paste0('--score /dcs04/nilanjan/data/jjin/prs/sim/results/ldpred2/ldpred2prs-mega-',race,'-rho=',rho,
                       '-size=',size,'-chr',chr,'-rep',rep,'-GA',GA,'.txt'),
                'cols=+scoresums,-scoreavgs',
                paste0('--score-col-nums 3-',ncol(beta_grid)),
                paste0(' --bfile /dcs04/nilanjan/data/jjin/prs/sim/',race,'/geno/mega/chr',chr),
                #paste0(' --bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/all_data/',race,'/chr',chr),
                " --threads 1",
                paste0(' --out /dcs04/nilanjan/data/jjin/prs/sim/results/ldpred2/ldpred2prs-mega-',race,'-rho=',rho,
                       '-size=',size,'-chr',chr,'-rep',rep,'-GA',GA))
system(prscode)


print(paste0('Complete training ldpred2 for ',race))
rm(beta_grid,corr0)



# ----------- Validation: run locally -----------
jindir = "/dcl01/chatterj/data/jin/"
RACES = c('EUR', 'AFR', 'AMR', 'EAS', 'SAS')
n.training = c(15000, 45000, 80000, 100000)
n.total = 1.2e5
split = list()
for (ethnicity in RACES){
  split[[ethnicity]] = list()
  for (siz in 1:4){
    split[[ethnicity]][[siz]] = list()
    n.test = n.val = (n.total-n.training[siz])/2#n.training[siz]*0.1
    sample.test = (n.training[siz]+1):(n.training[siz]+n.val)
    sample.val = (n.training[siz]+n.val+1):(n.total)
    split[[ethnicity]][[siz]] = list(1:n.training[siz], sample.test, sample.val)
  }
}
# ind.rows 
maindirec = "/dcl01/chatterj/data/yzhang/"
jindir = "/dcl01/chatterj/data/jin/"
output_path = paste0('/dcl01/chatterj/data/jin/prs/simulation/results/')
writescore_path = paste0('/dcs04/nilanjan/data/jjin/prs/sim/ldpred2/')

settings1 = rbind(expand.grid(rho = 1:3, size = 4, GA = 3:5, race = c('EUR')),
                  expand.grid(rho = 1:3, size = 1:4, GA = 3:5, race = c('AFR','AMR','EAS','SAS')))
settings2 = rbind(expand.grid(rho = 1:3, size = 4, GA = 1:2, race = c('EUR')),
                  expand.grid(rho = 1:3, size = 1:4, GA = 1:2, race = c('AFR','AMR','EAS','SAS')))
settings = rbind(settings1, settings2)
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))

R2.ldpred2 = cbind(settings,0,0,0,0,0)
colnames(p.tuned)[5:9] = c('R2','P-value','p','p.common','h2')

for (set in 1:nrow(settings)){
  race = as.character(settings[set,'race'])
  rho = settings[set,'rho']
  size = size2 = settings[set,'size']
  GA = settings[set,'GA']
  
  for (rep in 1:10){
    n.snp = numeric()
    #---------------------------------------#---------------------------------------
    chr = 22
    temfile = paste0('/dcs04/nilanjan/data/jjin/prs/sim/results/ldpred2/ldpred2prs-mega-',race,'-rho=',rho,
                     '-size=',size,'-chr',chr,'-rep',rep,'-GA',GA,'.txt')
    if(file.exists(temfile)){
      dftem1 = bigreadr::fread2(temfile)
      dftem1 = dftem1[,1:nrow(sets)]
      dftem1[is.na(dftem1)] = 0
      #dftem1 = dftem1 * (n.snp[chr]-1) * 2
      print(paste0('Chr ', chr,' Completed'))
    }
    if(!file.exists(temfile)) print(chr)
    for(chr in c(21:1)){
      temfile = paste0('/dcs04/nilanjan/data/jjin/prs/sim/results/ldpred2/ldpred2prs-mega-',race,'-rho=',rho,
                       '-size=',size,'-chr',chr,'-rep',rep,'-GA',GA,'.txt')
      if(file.exists(temfile)){
        a = bigreadr::fread2(temfile)
        a = a[,1:nrow(sets)]
        a[is.na(a)] = 0
        dftem1 = dftem1 + a #* n.snp[chr] * 2
        print(paste0('Chr ', chr,' Completed'))
      }
      if(!file.exists(temfile)) print(chr)
    }
    if (exists('dftem1')){
      write.table(dftem1,paste0(writescore_path,'mega-',race,'-rho=',rho,
                                '-size=',size,'-rep',rep,'-GA',GA,'.txt'), row.names = F,col.names = F, quote = FALSE, sep = "\t" )
      rm(dftem1)
    }
    
    validatetable = read.table(paste0('/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/',race,'/pheno_summary_out_GA/phenotypes_rho',rho,'_',GA,'.phen'),header=F,
                               colClasses = c(rep('character',2),rep('numeric',rep),rep('NULL',100-rep)))
    validatetable = as.data.frame(validatetable[,c(1,rep+2)])
    colnames(validatetable) = c('id','y')
    validatetable = validatetable[sort(unlist(split[[race]][[size]][2:3])),]
    
    output.train = matrix(0,nrow(sets),2); Output.validation = matrix(0,1,2)
    colnames(output.train) = colnames(Output.validation) = paste(as.vector(sapply(c('ldpred'),function(x){paste(c('R2','P'),x)})),race)
    rownames(output.train) = sapply(1:nrow(sets),function(x){paste0('p=',sets[x,1],' h2=', sets[x,2], 'sparse=', sets[x,3])})
    output.validation = Output.validation
    
    tem = bigreadr::fread2(paste0(writescore_path,'mega-',race,'-rho=',rho,'-size=',size,'-rep',rep,'-GA',GA,'.txt'))
    tem = sapply(tem[,1:ncol(tem)], as.numeric)
    # testing
    #prstable = data.frame(prs = tem[1:(nrow(tem)/2),], y = validatetable[1:(nrow(tem)/2),2])
    #formula.prs=formula(paste0(paste(paste0('y',rep), paste(c('ldpred'),collapse="+"), sep='~')))
    #fit = lm(formula.prs, data=prstable)
    output.train[,paste0('R2 ldpred ',race)] = as.numeric(sapply(1:nrow(sets),function(x){cor(tem[1:(nrow(tem)/2),x],validatetable[1:(nrow(tem)/2),2])^2}))
    #output.train[i,paste0('P ldpred ',race)] = summary(fit)$coefficients['ldpred','Pr(>|t|)']
    indx = which.max(output.train[,paste0('R2 ldpred ',race)])

    validatetable$ldpred = tem[,indx]
    prstable = validatetable[(nrow(tem)/2+1):(nrow(tem)),]
    prstable$ldpred = scale(prstable$ldpred,center=T,scale=T)
    formula.prs=formula(paste0(paste(paste0('y'), paste(c('ldpred'),collapse="+"), sep='~')))
    fit = lm(formula.prs, data=prstable)
    output.validation[1,paste0('R2 ldpred ',race)] = (coefficients(fit)['ldpred'])^2/var(na.omit(prstable[,paste0('y')]))
    output.validation[1,paste0('P ldpred ',race)] = summary(fit)$coefficients['ldpred','Pr(>|t|)']
    
    R2.ldpred2[set,c('R2','P-value','p','h2')] = R2.ldpred2[set,c('R2','P-value','p','h2')] + c(as.numeric(output.validation[1,paste0(c('R2 ldpred ', 'P ldpred '),race)]),signif(as.numeric(sets[indx,c('p','h2')]),2))#c(Output.validation/100,sets[indx,])
  }
  print(set)
}
R2.ldpred2 = R2.ldpred2/10
save(R2.ldpred2, file = paste0("/dcs04/nilanjan/data/jjin/prs/sim/ldpred2/R2.ldpred2.RData"))


