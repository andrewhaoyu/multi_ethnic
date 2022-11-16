#!/usr/bin/env bash
#$ -N pca
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu


/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/merge_1000G_ukbb/chrall \
--pca 20 allele-wts \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/plink/chrall

