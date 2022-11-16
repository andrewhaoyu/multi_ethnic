#!/usr/bin/env bash
#$ -N qc1000G
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=5000G
#$ -m e
#$ -t 1-22
#$ -M jzhan218@jhu.edu


/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \
--keep-allele-order \
--bfile /dcs04/nilanjan/data/jzhang2/1000G/GRCh37/original_files/chr${SGE_TASK_ID} \
--extract /dcs04/nilanjan/data/jzhang2/UKBB/genotype/allrsid/chr${SGE_TASK_ID}.biallelic.rsid /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_rsid/chr${SGE_TASK_ID}.biallelic.rsid \
--snps-only \
--maf 0.005 \
--geno 0.01 \
--indep-pairwise 50 10 0.1 \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_qc/chr${SGE_TASK_ID}_snps

mkdir /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_qc

for i in {21..22}
do

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \
--keep-allele-order \
--bfile /dcs04/nilanjan/data/jzhang2/1000G/GRCh37/original_files/chr${i} \
--extract /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_rsid/1000GUKBBoverlap.rsid \
--snps-only \
--maf 0.005 \
--geno 0.01 \
--indep-pairwise 50 10 0.1 \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_qc/chr${i}_snps

done


