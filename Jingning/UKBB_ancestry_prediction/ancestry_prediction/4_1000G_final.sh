#!/usr/bin/env bash
#$ -N final1000G
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=5000G
#$ -m e
#$ -M jzhan218@jhu.edu

mkdir /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_final/

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_qc/chrall \
--keep /dcs04/nilanjan/data/jzhang2/1000G/1000G_AFR_ID.txt /dcs04/nilanjan/data/jzhang2/1000G/1000G_AMR_ID.txt /dcs04/nilanjan/data/jzhang2/1000G/1000G_EAS_ID.txt /dcs04/nilanjan/data/jzhang2/1000G/1000G_EUR_ID.txt /dcs04/nilanjan/data/jzhang2/1000G/1000G_SAS_ID.txt \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_final/chrall

