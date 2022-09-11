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
#load 1kg snp list
setwd("/data/zhangh24/multi_ethnic/")
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
# input files
summary_EUR <- as.data.frame(fread(paste0(out.dir.sum,eth[1],"/summary_mega_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))


j = 1
summary_EUR_sub = summary_EUR %>% 
  filter(CHR==1) %>% 
  mutate(SE = BETA/Z) %>% 
  select(SNP, CHR, BP, effect_allele, non_effect_allele, BETA, SE, P, N)
#prepare the EUR data for SBayesR format
summary_EUR_XPASS = summary_EUR %>% 
  select(SNP, N, Z, effect_allele, non_effect_allele) %>% 
  rename(A1 = effect_allele,
         A2 = non_effect_allele)
write.table(summary_EUR_XPASS, file = paste0(temp.dir,eth[1],"summary_eur"),
            row.names = F, col.names = T, quote = T)