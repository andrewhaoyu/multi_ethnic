#!/usr/bin/env bash
#$ -N merge_final
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=5000G
#$ -m e
#$ -M jzhan218@jhu.edu

mkdir /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/merge_1000G_ukbb

echo /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_final/chrall >> /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/merge_1000G_ukbb/chrall.mergelist
echo /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_final/chrall >> /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/merge_1000G_ukbb/chrall.mergelist

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \
--merge-list /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/merge_1000G_ukbb/chrall.mergelist \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/merge_1000G_ukbb/chrall


