#!/usr/bin/env bash
#$ -N TG
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu


trait='TG'

## find the code of interested phenotype from ukbb html
#a <- read.table("/dcl01/arking/data/UK_Biobank/static/Phenotype/Downloads/tab_delim/ukb45830.tab", nrows=1)
#b <- as.character(a[1,])
#which(b == "f.30870.0.0")
#which(b == "f.30870.1.0")
mkdir /dcs04/nilanjan/data/jzhang2/UKBB/phenotype/${trait}
#cat /dcl01/arking/data/UK_Biobank/static/Phenotype/Downloads/tab_delim/ukb45830.tab | awk -F '\t' '{ print $1, $14610, $14611}' > /dcs04/nilanjan/data/jzhang2/UKBB/phenotype/${trait}/${trait}_from_ukbb.txt

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore "$1"  ../extract_phenotype_tuning+validation.R ${trait}.Rout
}
runr "--args trait='${trait}'"


