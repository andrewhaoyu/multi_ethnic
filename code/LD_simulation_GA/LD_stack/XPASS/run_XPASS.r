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

#####
# Install packages
if (!require("XPASS")) {
  if (!require("devtools")) {
    install.packages("devtools")
  }
  library(devtools)
  devtools::install_github("YangLabHKUST/XPASS")
}
packages <- c("data.table", "RhpcBLASctl", "getopt")
to_be_installed <- setdiff(packages, rownames(installed.packages()))
if (length(to_be_installed)) {
  install.packages(to_be_installed)
}

#####
# Load library
library(XPASS)
library(data.table)
library(RhpcBLASctl)
blas_set_num_threads(8)

#####
# parse the command line argument
eval(parse(text=paste(commandArgs(trailingOnly = TRUE), collapse=";")))

### if not using command line, then please specify the parameters below.
# race = "AFR" # target population
# size = 1 # size for target population
# GA = 1
# rho = 1
# rep = 1

#####
# input files
summary_EUR = paste0("/dcs04/nilanjan/data/wlu/XPASS/data/sumdata/EUR/sumdata-rho",rho,'-size4-rep',rep,'-GA',GA,'.txt') # auxilliary
summary_target = paste0("/dcs04/nilanjan/data/wlu/XPASS/data/sumdata/",race,'/sumdata-rho',rho,'-size',size,'-rep',rep,'-GA',GA,'.txt') # target
# auxilliary
ref_gene_EUR = "/dcs04/nilanjan/data/wlu/XPASS/data/ref_genotype/EUR/clump_ref_all_chr"
# target
ref_gene_target = paste0("/dcs04/nilanjan/data/wlu/XPASS/data/ref_genotype/", race, "/clump_ref_all_chr")
# predict
ref_gene_pred = paste0("/dcs04/nilanjan/data/wlu/XPASS/data/ref_genotype/", race, "/all_chr_test.mega")
# output file prefix
file_out = paste0("/dcs04/nilanjan/data/wlu/XPASS/XPASS_result/", race, "/XPASS-rho",rho,'-size',size,'-rep',rep,'-GA',GA)

#####
# XPASS
fit_bbj <- XPASS(file_z1 = summary_target, file_z2 = summary_EUR,
                 file_ref1 = ref_gene_target, file_ref2 = ref_gene_EUR,
                 # file_predGeno = ref_gene_pred, compPRS = T,
                 pop = "EUR", sd_method="LD_block", compPosMean = T,
                 file_out = file_out)

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





