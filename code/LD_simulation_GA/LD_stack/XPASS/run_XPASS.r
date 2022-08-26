# run_XPASS.R runs XPASS and calculate PRS scores (by both XPASS and plink,
# comment the last few lines out if plink is not needed). Note that the way
# to call plink might vary. (It's plink2 --... for me.)
#
# Run with:
# R CMD BATCH --vanill '--args race="AFR" size=1 GA=1 rho=2 rep=1' run_XPASS.R
#
# The section of "input files" controls the input and output file paths.
# Output files have the prefix paste0("/dcs04/nilanjan/data/wlu/XPASS/XPASS_result/", race, "/XPASS-rho",rho,'-size',size,'-rep',rep,'-GA',GA)

## clear workspace
rm(list = ls())
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])
i1 = 1
library(data.table)
library(tidyverse)
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/LD/'),showWarnings = FALSE)
temp.dir.LD <- paste0('/lscratch/',sid,'/test/LD/')
system(paste0("cp ",cur.dir,eth[1],"/clump_ref_all_chr.bed ",temp.dir,eth[1],"clump_ref_all_chr.bed"))
system(paste0("cp ",cur.dir,eth[1],"/clump_ref_all_chr.bim ",temp.dir,eth[1],"clump_ref_all_chr.bim"))
system(paste0("cp ",cur.dir,eth[1],"/clump_ref_all_chr.fam ",temp.dir,eth[1],"clump_ref_all_chr.fam"))
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.bed ",temp.dir,eth[i],"clump_ref_all_chr.bed"))
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.bim ",temp.dir,eth[i],"clump_ref_all_chr.bim"))
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.fam ",temp.dir,eth[i],"clump_ref_all_chr.fam"))

#####
# Load library
library(XPASS)
library(data.table)
library(RhpcBLASctl)
#blas_set_num_threads(8)

#####
# parse the command line argument
#eval(parse(text=paste(commandArgs(trailingOnly = TRUE), collapse=";")))

### if not using command line, then please specify the parameters below.
# race = "AFR" # target population
# size = 1 # size for target population
# GA = 1
# rho = 1
# rep = 1
i = 1
l = 1
m = 4
i_rep =1 
i1 = 1
#####
# input files
summary_EUR <- as.data.frame(fread(paste0(out.dir.sum,eth[1],"/summary_mega_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))
#summary_EUR = paste0("/dcs04/nilanjan/data/wlu/XPASS/data/sumdata/EUR/sumdata-rho",rho,'-size4-rep',rep,'-GA',GA,'.txt') # auxilliary
summary_target = as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_mega_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))) 
#summary_target = paste0("/dcs04/nilanjan/data/wlu/XPASS/data/sumdata/",race,'/sumdata-rho',rho,'-size',size,'-rep',rep,'-GA',GA,'.txt') # target
# auxilliary
ref_gene_EUR = paste0(temp.dir,eth[1],"clump_ref_all_chr")
# target
ref_gene_target = paste0(temp.dir,eth[i],"clump_ref_all_chr")
file_out = paste0("data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/xpass/rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)

#####
# XPASS
fit_bbj <- XPASS(file_z1 = summary_target, file_z2 = summary_EUR,
                 file_ref1 = ref_gene_target, file_ref2 = ref_gene_EUR,
                 # file_predGeno = ref_gene_pred, compPRS = T,
                 pop = "EUR", sd_method="LD_block", compPosMean = T,
                 file_out = file_out)

# predict
ref_gene_pred = paste0("/dcs04/nilanjan/data/wlu/XPASS/data/ref_genotype/", race, "/all_chr_test.mega")
# output file prefix




# separate in two steps due to huge REM needed
save(fit_bbj, file=paste0(file_out, "_param.RData"))

# heritability
heritability <- fit_bbj$H
write.table(heritability,file=paste0(file_out,"_Heritability.txt"),col.names = T,row.names = F,quote=F,sep="\t")

## Note:
# mu can be loaded afterward for PRS calculation using
# mu = fread(paste0(file_out, "_PosteriorMean.txt"))
## clear workspace
mu = fit_bbj$mu
rm(list = setdiff(ls(), c("mu", "ref_gene_pred", "file_out")))

# calculate PRS using XPASS
PRS <- predict_XPASS(mu,ref_gene_pred)
write.table(PRS,file=paste0(file_out,"_PRS.txt"),col.names = T,row.names = F,quote=F,sep="\t")

## clear workspace
rm(list = setdiff(ls(), c("mu", "ref_gene_pred", "file_out")))

# calculate PRS using plink (optional, comment the following out if not needed)
mu = mu[, c("SNP", "A1", "mu1", "mu2", "mu_XPASS1", "mu_XPASS2")]
write.table(mu,file = paste0(file_out,"_prs_prep"),col.names = T,row.names = F,quote=F)
res = system(paste0("plink2 --score-col-nums 3,4,5,6 --threads 2 --score ",file_out,"_prs_prep cols=+scoresums,-scoreavgs header no-mean-imputation --bfile ",ref_gene_pred," --out ",file_out,"_PRS_plink"))





