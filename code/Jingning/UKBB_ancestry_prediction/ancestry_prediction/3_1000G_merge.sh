#!/usr/bin/env bash
#$ -N merge1000G
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=5000G
#$ -m e
#$ -M jzhan218@jhu.edu


for i in {1..22}
do
echo /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_qc/chr${i} >> /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_qc/chrall.mergelist
done

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \
--keep-allele-order \
--merge-list /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_qc/chrall.mergelist \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_qc/chrall

