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
system(paste0("cp ",cur.dir,eth[1],"/clump_ref_all_chr.bed ",temp.dir,eth[1],"clump_ref_all_chr.bed"))
system(paste0("cp ",cur.dir,eth[1],"/clump_ref_all_chr.bim ",temp.dir,eth[1],"clump_ref_all_chr.bim"))
system(paste0("cp ",cur.dir,eth[1],"/clump_ref_all_chr.fam ",temp.dir,eth[1],"clump_ref_all_chr.fam"))
system(paste0("cp ",cur.dir,eth[1],"/clump_ref_all_chr_pca.eigenvec ",temp.dir,eth[1],"clump_ref_all_chr_pca.eigenvec"))
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.bed ",temp.dir,eth[i],"clump_ref_all_chr.bed"))
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.bim ",temp.dir,eth[i],"clump_ref_all_chr.bim"))
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr_pca.eigenvec ",temp.dir,eth[i],"clump_ref_all_chr_pca.eigenvec"))

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
#####
# input files
summary_EUR <- as.data.frame(fread(paste0(out.dir.sum,eth[1],"/summary_mega_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))
#prepare the EUR data for XPASS format
summary_EUR_XPASS = summary_EUR %>% 
  select(SNP, N, Z, effect_allele, non_effect_allele) %>% 
  rename(A1 = effect_allele,
         A2 = non_effect_allele)
write.table(summary_EUR_XPASS, file = paste0(temp.dir,eth[1],"summary_eur"),
            row.names = F, col.names = T, quote = F)
#prepare the pca of reference data for XPASS format
eur_ref_cov = as.data.frame(fread(paste0(temp.dir,eth[1],"clump_ref_all_chr_pca.eigenvec")))
eur_ref_cov = eur_ref_cov[,3:22]
write.table(eur_ref_cov, file = paste0(temp.dir,eth[1],"eur_ref_cov"),
            row.names = F, col.names = F, quote = F)

#summary_EUR = paste0("/dcs04/nilanjan/data/wlu/XPASS/data/sumdata/EUR/sumdata-rho",rho,'-size4-rep',rep,'-GA',GA,'.txt') # auxilliary
summary_target = as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_mega_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))) 
summary_target_XPASS = summary_target %>% 
  select(SNP, N, Z, effect_allele, non_effect_allele) %>% 
  rename(A1 = effect_allele,
         A2 = non_effect_allele)
write.table(summary_target_XPASS, file = paste0(temp.dir,eth[i],"summary_tar"),
            row.names = F, col.names = T, quote = F)
tar_ref_cov = as.data.frame(fread(paste0(temp.dir,eth[i],"clump_ref_all_chr_pca.eigenvec")))
tar_ref_cov = tar_ref_cov[,3:22]
write.table(tar_ref_cov, file = paste0(temp.dir,eth[i],"tar_ref_cov"),
            row.names = F, col.names = F, quote = F)

path_to_summary_eur = paste0(temp.dir,eth[1],"summary_eur")
path_to_summary_tar = paste0(temp.dir,eth[i],"summary_tar")
#summary_target = paste0("/dcs04/nilanjan/data/wlu/XPASS/data/sumdata/",race,'/sumdata-rho',rho,'-size',size,'-rep',rep,'-GA',GA,'.txt') # target
# auxilliary
ref_gene_EUR = paste0(temp.dir,eth[1],"clump_ref_all_chr")
# target
ref_gene_target = paste0(temp.dir,eth[i],"clump_ref_all_chr")

path_to_cov_EUR = paste0(temp.dir,eth[1],"eur_ref_cov")
path_to_cov_tar = paste0(temp.dir,eth[i],"tar_ref_cov")

file_out = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/xpass/cov_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)

#####
# XPASS
fit_bbj <- XPASS(file_z1 = path_to_summary_tar, file_z2 = path_to_summary_eur,
                 file_ref1 = ref_gene_target, file_ref2 = ref_gene_EUR,
                 file_cov1 = path_to_cov_tar,file_cov2 = path_to_cov_EUR,
                 # file_predGeno = ref_gene_pred, compPRS = T,
                 pop = "EUR", sd_method="LD_block", compPosMean = T,
                 file_out = file_out)
save(fit_bbj, file=paste0(file_out, "_param.RData"))
print("XPASS finished run")

system(paste0("cp ",cur.dir,eth[i],"/all_chr_test.mega.* ",temp.dir))

# predict
ref_gene_pred = paste0(temp.dir,"/all_chr_test.mega")
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





