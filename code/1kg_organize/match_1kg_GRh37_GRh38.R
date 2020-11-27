#match 1kg GRCh37 and GRCh38
#download vcf file for 1kg https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/
#delete teh first 56 document lines sed -e '1,56d' 00-common_all.vcf > 1kg_clean.vcf
#keep the rsid and chr:pos: " awk '{print $3, $1":"$2}' 1kg_clean.vcf > 1kg_clean_update.vcf"
#only keep 1kg_clean_update.vcf "mv 1kg_clean_update.vcf 1kg_clean.vcf"
library(bigsnpr)
library(bigreadr)
library(data.table)
library(dplyr)
vcf.infor <- as.data.frame(fread2("/gpfs/gsfs11/users/zhangh24/KG.vcf/1kg_clean.vcf",header=T))
colnames(vcf.infor) <- c("rs_id","chr.pos")
#colnames(vcf.infor)[1] <- "chr"
#idx <- which(vcf.infor$POS==10642)
setwd("/data/zhangh24/multi_ethnic/")
#load 1kg snp list
load("./result/LD_simulation_new/snp.infor.rdata")
snp.infor.new = snp.infor %>% 
  mutate(chr.pos = paste0(CHR,":",position))
snp.infor.update = left_join(snp.infor.new,vcf.infor,by="chr.pos")

#save(snp.infor.update,file = "./result/LD_simulation_new/snp.infor.update.rdata")
load("./result/LD_simulation_new/snp.infor.38.rdata")

snp.infor.match = left_join(snp.infor.update,
                            snp.infor.38,by="rs_id") %>% 
  mutate(position_GRCh37=position) 

idx <- which(snp.infor.match$position_GRCh37!=snp.infor.match$position_GRCh38)

snp.infor.match = snp.infor.match %>% 
  select(id,position,a0,a1,TYPE,AFR,AMR,EAS,EUR,SAS,ALL,CHR,rs_id,position_GRCh38,position_GRCh37)
save(snp.infor.match,file = "./result/LD_simulation_new/snp.infor.match37_38.rdata")
