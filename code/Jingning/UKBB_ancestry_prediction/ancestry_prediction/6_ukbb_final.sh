#!/usr/bin/env bash
#$ -N ukbb_final
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=5000G
#$ -m e
#$ -M jzhan218@jhu.edu


for i in {1..22}
do
echo /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_qc/chr${i} >> /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_qc/chrall.mergelist
done

mkdir /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_final
/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \
--keep-allele-order \
--merge-list /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_qc/chrall.mergelist \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_final/chrall


