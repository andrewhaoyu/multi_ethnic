#use PRS-csx on simulation
library(data.table)
library(dplyr)
setwd("/data/zhangh24/test1/ref1/")
snp_mult = fread("snpinfo_mult_1kg_hm3")

snp_mult = snp_mult %>% 
  filter(CHR==22)
write.table(snp_mult, file = "./snpinfo_mult_1kg_hm3",row.names = F, col.names = T, quote = F)
idx <- which(snp_mult$SNP=="rs174342")
snp_mult[idx,]
snpinfor  = fread("./ldblk_1kg_eur/snpinfo_1kg_hm3")
snpinfor = snpinfor %>% 
  filter(CHR == 22)
write.table(snpinfor, file = "./ldblk_1kg_eur/snpinfo_1kg_hm3",row.names = F, col.names = T, quote = F)
idx <- which(snpinfor$SNP=="rs174342")
snpinfor[idx,]
snpinfor  = fread("./ldblk_1kg_eas/snpinfo_1kg_hm3")
snpinfor = snpinfor %>% 
  filter(CHR == 22)
write.table(snpinfor, file = "./ldblk_1kg_eas/snpinfo_1kg_hm3",row.names = F, col.names = T, quote = F)
idx <- which(snpinfor$SNP=="rs174342")
print(idx)
snpinfor[idx,]






setwd("/data/zhangh24/test1/ref2/")
snp_mult = fread("snpinfo_mult_1kg_hm3")

snp_mult = snp_mult %>% 
  filter(CHR==22)
write.table(snp_mult, file = "./snpinfo_mult_1kg_hm3",row.names = F, col.names = T, quote = F, sep= "\t")
idx <- which(snp_mult$SNP=="rs174342")
snp_mult[idx,]
snpinfor  = fread("./ldblk_1kg_eur/snpinfo_1kg_hm3")
snpinfor = snpinfor %>% 
  filter(CHR == 22)
write.table(snpinfor, file = "./ldblk_1kg_eur/snpinfo_1kg_hm3",row.names = F, col.names = T, quote = F, sep= "\t")
idx <- which(snpinfor$SNP=="rs174342")
snpinfor[idx,]
snpinfor  = fread("./ldblk_1kg_eas/snpinfo_1kg_hm3")
snpinfor = snpinfor %>% 
  filter(CHR == 22)
write.table(snpinfor, file = "./ldblk_1kg_eas/snpinfo_1kg_hm3",row.names = F, col.names = T, quote = F, sep = "\t")

idx <- which(snpinfor$SNP=="rs174342")
snpinfor[idx,]


