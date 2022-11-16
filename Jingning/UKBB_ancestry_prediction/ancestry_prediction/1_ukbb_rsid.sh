#!/usr/bin/env bash
#$ -N rsid
#$ -cwd
#$ -l mem_free=1G,h_vmem=1G,h_fsize=100G
#$ -m e
#$ -t 1-22
#$ -M jzhan218@jhu.edu


cat /dcl01/arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr${SGE_TASK_ID}.pvar | awk -F '\t' 'NR>1{print $3}' > /dcs04/nilanjan/data/jzhang2/UKBB/genotype/allrsid/chr${SGE_TASK_ID}.rsid

mkdir /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_rsid/

R

gnomad <- readLines("/dcs04/nilanjan/data/jzhang2/UKBB/*backup_gnomad_v3/1_subset_gnomad_variants/gnomad_variats_v3_rsid.txt")

for(chr in 1:22){
  print(chr)
  a <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/genotype/allrsid/chr",chr,".rsid"))
  a <- a[a %in% gnomad]
  b <- table(a)
  c <- names(b[b>1])
  writeLines(a[!(a %in% c)], paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_rsid/chr",chr,".biallelic.rsid"))
}

res <- character()
for(chr in 1:22){
  print(chr)
  a1 <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/ukbb_rsid/chr",chr,".biallelic.rsid"))
  a2 <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_rsid/chr",chr,".biallelic.rsid"))
  b <- intersect(a1, a2)
  res <- c(res, b)
}
writeLines(res, paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_rsid/1000GUKBBoverlap.rsid"))





#a <- readLines("/dcs04/nilanjan/data/jzhang2/UKBB/genotype/allrsid/chrall.rsid")
#b <- table(a)
#names(b>1)
#sum(a=="rs1015011")

