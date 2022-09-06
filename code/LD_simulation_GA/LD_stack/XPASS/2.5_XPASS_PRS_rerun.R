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
i1 =as.numeric(args[[5]])
library(data.table)
library(tidyverse)
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')


# output file prefix
load(paste0(file_out, "_param.RData"))
mu = fit_bbj$mu
# calculate PRS using plink (optional, comment the following out if not needed)
mu = mu[, c("SNP", "A1", "mu1", "mu2", "mu_XPASS1", "mu_XPASS2")]
write.table(mu,file = paste0(temp.dir,"prs_prep"),col.names = T,row.names = F,quote=F)
res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                    "--score-col-nums 3,4,5,6 --threads 2 ",
                    "--score ",temp.dir,"prs_prep cols=+scoresums,-scoreavgs header no-mean-imputation ",
                    "--bfile ",ref_gene_pred,
                    " --out ",file_out,"_PRS"))





