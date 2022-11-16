#!/usr/bin/env bash
#$ -N ukbb_qc
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=5000G
#$ -t 1-22
#$ -m e
#$ -M jzhan218@jhu.edu

mkdir /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_qc

#R

#a <- read.table("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/unrelated_whites.id",header=T)
#set.seed(1)
#b <- a[sample(1:nrow(a), 8000),]
#write_tsv(b, "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_qc/unrelated_whites_8000.id")
#
#a <- read.table("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/self_reported_White.id", header=T)
#b <- a[!(a$IID %in% b$IID),]
#write_tsv(b, "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_qc/self_reported_White_exclude_8000_unrelated_whites.id")

#cat /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_final/chrall.bim | awk -F '\t' '{print $2}' > /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_final/chrall.rsid


#for i in {1..22}
#do
/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--pfile /dcl01/arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr$SGE_TASK_ID \
--extract /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_final/chrall.rsid \
--remove /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_qc/self_reported_White_exclude_8000_unrelated_whites.id \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_qc/chr$SGE_TASK_ID
#done

