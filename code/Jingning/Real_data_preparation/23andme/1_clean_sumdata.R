

disease='any_cvd'
race='EUR'


library(readr)
library(bigreadr)
library(dplyr)

## get snp info of mega snps (mega)

gt <- fread2(paste0("/dcl01/chatterj/data/23andme/",race,"/annotations/",race,"/gt_snp_stat.txt"))
im <- fread2(paste0("/dcl01/chatterj/data/23andme/",race,"/annotations/",race,"/im_snp_stat.txt"))
mega <- readLines("/dcl01/chatterj/data/jin/MEGA/megarsid.txt")
gt <- gt[gt$assay.name %in% mega,]
write_tsv(gt, paste0("/dcl01/chatterj/data/23andme/",race,"/annotations/gt_snp_stat_mega.txt"))
im <- im[im$assay.name %in% mega,]
write_tsv(im, paste0("/dcl01/chatterj/data/23andme/",race,"/annotations/im_snp_stat_mega.txt"))


########################################################################
########################################################################
# clean summary data (mega)

dir.create(paste0("/dcl01/chatterj/data/23andme/",race,"/sumdat/cleaned/"))

mega <- readLines("/dcl01/chatterj/data/jin/MEGA/megarsid.txt")

a <- fread2(paste0("/dcl01/chatterj/data/23andme/",race,"/sumdat/",disease,".dat"))
b <- a[(!is.na(a$pvalue)) & (a$pass=="Y"),] # non-missing p-value and passing QC

# 1. snp in summary data:
dim(a)
# 2. snp in summary data passing QC and non-missing p-value:
dim(b)

write_tsv(b,paste0("/dcl01/chatterj/data/23andme/",race,"/sumdat/cleaned/",disease,"_passQC_noNA.dat"))

b0 <- b[,c(1:5,7:10)]
write_tsv(b0,paste0("/dcl01/chatterj/data/23andme/",race,"/sumdat/cleaned/",disease,"_passQC_noNA_23andmeid.dat"))

gt <- read_tsv(paste0("/dcl01/chatterj/data/23andme/",race,"/annotations/gt_snp_stat_mega.txt"))
im <- read_tsv(paste0("/dcl01/chatterj/data/23andme/",race,"/annotations/im_snp_stat_mega.txt"))

b1 <- b0[b0$src=="G",] ## get genotyped variants
b2 <- b0[b0$src=="I",] ## get imputed variants
b1 <- inner_join(b1, gt[,c(1,2,5)], by=c("all.data.id"="gt.data.id")) # match rsid for genotyped variants
b2 <- inner_join(b2, im[,c(1,2,4)], by=c("all.data.id"="im.data.id")) # match rsid for imputed variants
b_mega <- rbind(b1, b2)

write_tsv(b_mega,paste0("/dcl01/chatterj/data/23andme/",race,"/sumdat/cleaned/",disease,"_passQC_noNA_rsid_mega.dat"))

# 3. mega snp in summary data passing QC and non-missing p-value:
dim(b_mega)

# 4. mega snp in summary data passing QC and non-missing p-value with allele frequency:
sum(!is.na(b_mega$freq.b))

# 5. common mega snp in summary data passing QC and non-missing p-value with allele frequency:
sum(b_mega$freq.b<0.99 & b_mega$freq.b>0.01, na.rm=T)

mean(b_mega$freq.b<0.99 & b_mega$freq.b>0.01, na.rm=T)



library(bigreadr)
library(readr)
library(dplyr)

nsnps <- integer()
nsnps_mega <- integer()
i=0

for (disease in c("any_cvd","depression","heart_metabolic_disease_burden","height","iqb.sing_back_musical_note","migraine_diagnosis","morning_person")){
  for (race in c("EUR","AFR","AMR","EAS","SAS")){
    print(paste0(disease,"-",race))
    i=i+1
    if("RST-8714-im" %in% dir(paste0("/dcs04/nilanjan/data/23andme/",race,"/sumdat/"))){
      a <- fread2(paste0("/dcs04/nilanjan/data/23andme/",race,"/sumdat/RST-8714-im/",disease,".dat"))
    }else{
      a <- fread2(paste0("/dcs04/nilanjan/data/23andme/",race,"/sumdat/",disease,".dat"))
    }
    if("RST-8714-gt" %in% dir(paste0("/dcs04/nilanjan/data/23andme/",race,"/sumdat/"))){
      a2 <- fread2(paste0("/dcs04/nilanjan/data/23andme/",race,"/sumdat/RST-8714-gt/",disease,".dat"))
    }else{
      a2 <- fread2(paste0("/dcs04/nilanjan/data/23andme/",race,"/sumdat/",disease,".dat"))
    }


    b <- a[(!is.na(a$pvalue)),] # non-missing p-value
    nsnps[i] <- nrow(b)

    b0 <- b[,c(1:5,7:10)]
    gt <- read_tsv(paste0("/dcs04/nilanjan/data/23andme/annotations/",race,"_gt_snp_stat_mega.txt"))
    im <- read_tsv(paste0("/dcs04/nilanjan/data/23andme/annotations/",race,"_im_snp_stat_mega.txt"))
    b1 <- b0[b0$src=="G",] ## get genotyped variants
    b2 <- b0[b0$src=="I",] ## get imputed variants
    b1 <- inner_join(b1, gt[,c(1,2,5)], by=c("all.data.id"="gt.data.id")) # match rsid for genotyped variants
    b2 <- inner_join(b2, im[,c(1,2,4)], by=c("all.data.id"="im.data.id")) # match rsid for imputed variants
    b_mega <- rbind(b1, b2)
    nsnps_mega[i] <- nrow(b_mega)

  }
}
res <- data.frame(nsnps, nsnps_mega)
res



