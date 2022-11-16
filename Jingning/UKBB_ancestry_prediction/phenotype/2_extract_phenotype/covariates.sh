#!/usr/bin/env bash
#$ -N covariates
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu


trait='covariates'

## find the code of interested phenotype from ukbb html
#a <- read.table("/dcl01/arking/data/UK_Biobank/static/Phenotype/Downloads/tab_delim/ukb45830.tab", nrows=1)
#b <- as.character(a[1,])
#which(b == "f.31.0.0") # sex: 23
#which(b == "f.21003.0.0") # age: 9773
#which(b == "f.22009.0.1") # PC1: 9948
#which(b == "f.22009.0.10") # PC10: 9957

mkdir /dcs04/nilanjan/data/jzhang2/UKBB/phenotype/${trait}
#cat /dcl01/arking/data/UK_Biobank/static/Phenotype/Downloads/tab_delim/ukb45830.tab | awk -F '\t' '{ print $1, $23, $9773, $9948, $9949, $9950, $9951, $9952, $9953, $9954, $9955, $9956, $9957}' > /dcs04/nilanjan/data/jzhang2/UKBB/phenotype/${trait}/${trait}_from_ukbb.txt

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore "$1"  ../extract_phenotype_tuning+validation.R ${trait}.Rout
}
runr "--args trait='${trait}'"





