#!/usr/bin/env bash
#$ -N ukbb_qc
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=5000G
#$ -m e
#$ -M jzhan218@jhu.edu

mkdir /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_qc

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--pfile /dcl01/arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr21 \
--extract /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/pc1000G_5.snps \
--remove /dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/self_reported_White.id \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_qc/chr21

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--pfile /dcl01/arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr22 \
--extract /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/pc1000G_5.snps \
--remove /dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/self_reported_White.id \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_qc/chr22

for i in {21..22}
do
echo /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_qc/chr${i} >> /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_qc/chrall.mergelist
done

mkdir /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_final
/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \
--keep-allele-order \
--merge-list /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_qc/chrall.mergelist \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_final/chr2122


