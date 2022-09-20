#goal: calculate the top 20PCs for the 1KG data
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
library(data.table)
library(tidyverse)
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
kg.dir = "/data/zhangh24/KGref_MEGA/GRCh37/"
system(paste0("cp ",kg.dir,eth[i],"/all_chr.bed ",temp.dir,eth[i],"all_chr.bed"))
system(paste0("cp ",kg.dir,eth[i],"/all_chr.bim ",temp.dir,eth[i],"all_chr.bim"))
system(paste0("cp ",kg.dir,eth[i],"/all_chr.fam ",temp.dir,eth[i],"all_chr.fam"))


ref_gene = paste0(temp.dir,eth[i],"all_chr")
file_out = paste0(temp.dir,"all_chr_pca")
res = system(paste0("/data/zhangh24/software/plink2 ",
                    "--threads 2 ",
                    "--pca ",
                    "--bfile ",ref_gene,
                    " --out ",file_out))
system(paste0("mv ",file_out,".eigenvec ",
              kg.dir,eth[i],"/"))
