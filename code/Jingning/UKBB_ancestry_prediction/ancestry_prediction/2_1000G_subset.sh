#!/usr/bin/env bash
#$ -N subset1000G
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=5000G
#$ -m e
#$ -t 1-22
#$ -M jzhan218@jhu.edu

mkdir /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_qc
for i in {1..22}
do
/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--bfile /dcs04/nilanjan/data/jzhang2/1000G/GRCh37/original_files/chr${i} \
--extract /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_rsid/1000GUKBBoverlap.rsid \
--rm-dup exclude-all \
--snps-only \
--maf 0.001 \
--geno 0.01 \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_qc/chr${i}
done

