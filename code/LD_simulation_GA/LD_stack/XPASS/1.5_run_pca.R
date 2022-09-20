#goal: calculate the top 20PCs for the 
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
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.bed ",temp.dir,eth[i],"clump_ref_all_chr.bed"))
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.bim ",temp.dir,eth[i],"clump_ref_all_chr.bim"))
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.fam ",temp.dir,eth[i],"clump_ref_all_chr.fam"))
ref_gene = paste0(temp.dir,eth[i],"clump_ref_all_chr")
file_out = paste0(temp.dir,"clump_ref_all_chr_pca")
res = system(paste0("/data/zhangh24/software/plink2 ",
                    "--threads 2 ",
                    "--pca ",
                    "--bfile ",ref_gene,
                    " --out ",file_out))
system(paste0("mv ",file_out,".eigenvec ",
              cur.dir,eth[i],"/"))
#temp = fread("/lscratch/48000777/test/clump_ref_all_chr_pca.eigenvec")

