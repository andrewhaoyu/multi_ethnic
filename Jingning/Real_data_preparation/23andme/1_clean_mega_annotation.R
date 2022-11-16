
library(readr)
library(bigreadr)
library(dplyr)

race='AMR'

## get snp info of mega snps (mega)

mega <- readLines("/dcl01/chatterj/data/jin/MEGA/megarsid.txt")

gt <- fread2(paste0("/dcl01/chatterj/data/23andme/",race,"/annotations/",race,"/gt_snp_stat.txt"))
save(gt, file = paste0("/dcl01/chatterj/data/23andme/",race,"/annotations/gt_snp_stat.RData"))
gt <- gt[gt$assay.name %in% mega,]
write_tsv(gt, paste0("/dcl01/chatterj/data/23andme/",race,"/annotations/gt_snp_stat_mega.txt"))
rm(list = "gt")
im <- fread2(paste0("/dcl01/chatterj/data/23andme/",race,"/annotations/",race,"/im_snp_stat.txt"))
save(im, file = paste0("/dcl01/chatterj/data/23andme/",race,"/annotations/im_snp_stat.RData"))
im <- im[im$assay.name %in% mega,]
write_tsv(im, paste0("/dcl01/chatterj/data/23andme/",race,"/annotations/im_snp_stat_mega.txt"))
rm(list = "im")

snpinfo <- fread2("/dcl01/chatterj/data/23andme/snpinfo/all_snp_info.txt")
save(snpinfo, file = paste0("/dcl01/chatterj/data/23andme/snpinfo/all_snp_info.RData"))




mega <- readLines("/dcl01/chatterj/data/jin/MEGA/mega-hm3-rsid.txt")

for (race in c("EUR","AFR","AMR","EAS","SAS")){
  gt <- fread2(paste0("/dcs04/nilanjan/data/23andme/",race,"/annotations/",race,"/gt_snp_stat.txt"))
  #save(gt, file = paste0("/dcs04/nilanjan/data/23andme/",race,"/annotations/gt_snp_stat.RData"))
  gt <- gt[gt$assay.name %in% mega,]
  write_tsv(gt, paste0("/dcs04/nilanjan/data/23andme/annotations/",race,"_gt_snp_stat_mega.txt"))
  rm(list = "gt")
  im <- fread2(paste0("/dcs04/nilanjan/data/23andme/",race,"/annotations/",race,"/im_snp_stat.txt"))
  #save(im, file = paste0("/dcs04/nilanjan/data/23andme/",race,"/annotations/im_snp_stat.RData"))
  im <- im[im$assay.name %in% mega,]
  write_tsv(im, paste0("/dcs04/nilanjan/data/23andme/annotations/",race,"_im_snp_stat_mega.txt"))
  rm(list = "im")
}
#snpinfo <- fread2("/dcs04/nilanjan/data/23andme/snpinfo/all_snp_info.txt")
#save(snpinfo, file = paste0("/dcs04/nilanjan/data/23andme/snpinfo/all_snp_info.RData"))

