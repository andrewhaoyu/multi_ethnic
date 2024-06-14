rm(list = ls())

library(dplyr)
library(stringr)
library(readxl)


## /data/BB_Bioinformatics/ProjectData/1000G_full_data/GRCh37/


eth <- c("EUR","AFR","EAS","AMR","SAS")

array <- as.numeric(commandArgs(TRUE)[1])
if(array < 23){
  i <- 1
  j <- array - 0
}else if(array < 45){
  i <- 2
  j <- array - 22
}else if(array < 67){
  i <- 3
  j <- array - 44
}else if(array < 89){
  i <- 4
  j <- array - 66
}else{
  i <- 5
  j <- array - 88
}

local_dir <- paste0("/data/williamsjacr/ldsc_corrected/",eth[i],"/")

if(i!=2){
  #AMR 1000 genomes data is relatively limited
  system(paste0("ldsc.py ",
                "--bfile ",local_dir,"chr_clean",j," ",
                "--print-snps /data/williamsjacr/ldsc_corrected/w_hm3_ldsc_printsnps.snplist"," ",
                "--l2 --ld-wind-kb 1000 ",
                "--out ",local_dir,eth[i],"_hm3_corrected_ldsc/",j))
}else{
  system(paste0("ldsc.py ",
                "--bfile ",local_dir,"chr_clean",j," ",
                "--print-snps /data/williamsjacr/ldsc_corrected/w_hm3_ldsc_printsnps.snplist"," ",
                "--l2 --ld-wind-kb 2000 ",
                "--out ",local_dir,eth[i],"_hm3_corrected_ldsc/",j))
}