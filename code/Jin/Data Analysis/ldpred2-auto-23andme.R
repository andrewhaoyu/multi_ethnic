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
traits = c('any_cvd','heart_metabolic_disease_burden','depression','height',
           'morning_person','migraine_diagnosis','iqb.sing_back_musical_note')
races = c('EUR','AFR','AMR','SAS','EAS')
temp <- commandArgs(TRUE)
race = races[as.numeric(temp[1])]
chr =  as.numeric(temp[2])
ldr = 1/1000
ncores = 1

#### read in reference data
temfile = paste0('/dcs04/nilanjan/data/jjin/prs/1KGref_MEGA/GRCh37/',race,'/chr',chr,'.bk')
system(paste0('rm -rf ',temfile))
temfile = paste0('/dcs04/nilanjan/data/jjin/prs/1KGref_MEGA/GRCh37/',race,'/chr',chr,'.rds')
system(paste0('rm -rf ',temfile))
#if (!file.exists(temfile)){
snp_readBed(paste0('/dcs04/nilanjan/data/jjin/prs/1KGref_MEGA/GRCh37/',race,'/chr',chr,'.bed'))
#}
obj.bigSNP <- snp_attach(paste0('/dcs04/nilanjan/data/jjin/prs/1KGref_MEGA/GRCh37/',race,'/chr',chr,'.rds'))
map <- obj.bigSNP$map[-c(3)]
names(map) <- c("chr", "rsid", "pos", "a0", "a1") # c("chr", "pos", "a0", "a1")

G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
POS2 <- snp_asGeneticPos(CHR, POS, dir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/23andme/',race,'/geno/'), ncores = ncores)
NCORES <-  8 #nb_cores()
#trait = 'height'
for (trait in traits){
  temdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/23andme/',race,'/score/')
  temfile = paste0(temdir,"ldpred-",trait,"-chr",chr,"-auto-ldr",ldr,".txt")
  temRData = paste0(temdir,"ldpred-",trait,"-chr",chr,"-auto-pars-ldr",ldr,".RData")
  if (file.exists(temRData)){
    print(paste0('Completed data preparation: ',trait))
    print(paste0('Completed: ',trait))
  }
  if (!file.exists(temRData)){
    # Read external summary statistics
    sumraw = bigreadr::fread2(paste0('/dcs04/nilanjan/data/23andme/cleaned/',race,'/sumdat/',trait,'_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt'))
    if (chr == 6){
      sumraw = sumraw[((sumraw$BP<=26e6) | (sumraw$BP>=34e6)),]
    }
    sumraw$n_eff = sumraw$N_control
    if (trait %in% c('any_cvd','depression','morning_person','migraine_diagnosis','iqb.sing_back_musical_note')){
      sumraw$n_eff = 4 / (1 / sumraw$N_case + 1 / sumraw$N_control)
    }
    sumraw$N_case <- sumraw$N_control <- NULL
    
    sumstats = sumraw[sumraw$CHR == chr,c('CHR','rsid','BP','A1','A2','BETA','SD','P','n_eff')]

    set.seed(2020)
    names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
    sumstats = sumstats[sumstats$rsid %in% map$rsid,]
    info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
    rownames(info_snp) = info_snp$rsid
    ## compute correlation
    ind.chr <- which(info_snp$chr == chr)
    df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
    ## indices in G
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    
    corr0 <- snp_cor(G, ind.col = ind.chr2, ncores = NCORES,
                     infos.pos = POS2[ind.chr2], size =  ldr)
    corr <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))
    
    # Automatic model
    ldsc <- snp_ldsc2(corr0, df_beta)
    h2_est <- abs(ldsc[["h2"]])
    print(paste0('Complete data preparation: ',trait))
    
    multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                                   vec_p_init = seq_log(1e-4, 0.5, length.out = 15),
                                   ncores = NCORES)
    save(multi_auto, file = temRData)
    beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
    beta_auto0 = cbind(beta_auto, info_snp[,c('a0','a1','rsid')])
    write_delim(beta_auto0,file = temfile, delim='\t')
    rm(corr0, corr)
    print(paste0('Complete: ',trait))
  }
}
