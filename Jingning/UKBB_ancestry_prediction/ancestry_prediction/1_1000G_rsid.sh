#!/usr/bin/env bash
#$ -N rsid
#$ -cwd
#$ -l mem_free=1G,h_vmem=1G,h_fsize=100G
#$ -m e
#$ -t 1-22
#$ -M jzhan218@jhu.edu

mkdir /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_rsid
for i in {1..22}
do
cat /dcs04/nilanjan/data/jzhang2/1000G/GRCh37/original_files/chr${i}.bim | awk -F '\t' '{print $2}' > /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_rsid/chr${i}.rsid
done

R

gnomad <- readLines("/dcs04/nilanjan/data/jzhang2/UKBB/*backup_gnomad_v3/1_subset_gnomad_variants/gnomad_variats_v3_rsid.txt")
for(chr in 1:22){
  print(chr)
  a <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_rsid/chr",chr,".rsid"))
  a <- a[a %in% gnomad]
  b <- table(a)
  c <- names(b[b>1])
  writeLines(a[!(a %in% c)], paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_rsid/chr",chr,".biallelic.rsid"))
}

#a <- readLines("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_rsid/chrall.rsid")
#b <- table(a)
#names(b>1)
#sum(a=="rs1015011")

