
#trait='HDL'
#race='AMR'

allele.qc = function(a1,a2,ref1,ref2) {
  a1 = toupper(a1)
  a2 = toupper(a2)
  ref1 = toupper(ref1)
  ref2 = toupper(ref2)

  ref = ref1
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip1 = flip

  ref = ref2
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip2 = flip;
  #
  #snp = list()
  #snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  #snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
  #snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
  #snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
  #snp[["amb_str"]] = (a1 == flip2 & a2 == flip1) | (a1 == flip1 & a2 == flip2)

  snp = list()
  snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
  snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
  snp[["keep"]][ !( (a1 == flip2 & a2 == flip1) | (a1 == flip1 & a2 == flip2) | (a1 == ref2 & a2 == ref1) | (a1 == ref1 & a2 == ref2) ) ] = F
  snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
  snp[["amb_str"]] = paste0(a1,":",a2) %in% c("A:T","T:A","G:C","C:G")
    #(a1 == flip2 & a2 == flip1) | (a1 == flip1 & a2 == flip2)

  return(snp)
}

library(readr)
library(bigreadr)
library(dplyr)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

num_rec <- integer()

# clean summary data (mega)

dir.create("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned")
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",race))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",race,"/clean_step/"))

mega <- readLines("/dcl01/chatterj/data/jin/MEGA/mega-hm3-rsid.txt")

# 1. variants in summary data:
a <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/GLGC_raw/",race,"/",trait))
num_rec[1] <- dim(a)[1]

# 2. mega snp
b <- a[a$rsID %in% mega,]
num_rec[2] <- dim(b)[1]; rm("a")

# 3. snp in summary data non-missing p-value:
b_mega <- b[(!is.na(b$pvalue)) & (nchar(b$REF)==1) & (nchar(b$ALT)==1),] # non-missing p-value SNP
num_rec[3] <- dim(b_mega)[1]; rm("b")

# 4. remove duplicate SNP
tmp <- table(b_mega$rsID)
dup <- names(tmp[tmp>1])
b_mega_nodup <- b_mega[! b_mega$rsID %in% dup, ]
num_rec[4] <- dim(b_mega_nodup)[1]; rm("b_mega")
write_tsv(b_mega_nodup, paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",race,"/clean_step/",trait,"_noNA_SNP_mega_nodup.dat"))

# 5. af
af <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/1000G/GRCh37/",race,"/mega_freq.afreq"))
b_mega_nodup_af <- left_join(b_mega_nodup, af[,c(2:5)], by=c("rsID"="ID"))
## in 1000G
b_mega_nodup_af <- b_mega_nodup_af[!is.na(b_mega_nodup_af$ALT_FREQS),]
num_rec[5] <- dim(b_mega_nodup_af)[1]; rm("b_mega_nodup")
qc <- allele.qc(a1=b_mega_nodup_af$ALT.x, a2=b_mega_nodup_af$REF.x, ref1=b_mega_nodup_af$ALT.y, ref2=b_mega_nodup_af$REF.y)

# 6. remove ambiguous strand
b_mega_nodup_af$A1 <- b_mega_nodup_af$ALT.y
b_mega_nodup_af$A2 <- b_mega_nodup_af$REF.y
b_mega_nodup_af$A1_AF_1000G <- b_mega_nodup_af$ALT_FREQS
b_mega_nodup_af$POOLED_ALT_AF[qc$flip] <- 1 - b_mega_nodup_af$POOLED_ALT_AF[qc$flip]
b_mega_nodup_af$EFFECT_SIZE[qc$flip] <- -b_mega_nodup_af$EFFECT_SIZE[qc$flip]
b_mega_nodup_af_noamb <- b_mega_nodup_af[qc$keep,]
num_rec[6] <- dim(b_mega_nodup_af_noamb)[1]; rm("b_mega_nodup_af")

# 7. remove small sample size
N = median(b_mega_nodup_af_noamb$N)
b_mega_nodup_af_noamb_rmN <- b_mega_nodup_af_noamb[b_mega_nodup_af_noamb$N > N*0.9, ]
num_rec[7] <- dim(b_mega_nodup_af_noamb_rmN)[1]; rm("b_mega_nodup_af_noamb")
write_tsv(b_mega_nodup_af_noamb_rmN, paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",race,"/clean_step/",trait,"_noNA_SNP_mega_nodup_noambiguous_nosmallN.dat"))

# 8. common
b_mega_nodup_af_noamb_rmN_comm <- b_mega_nodup_af_noamb_rmN[(b_mega_nodup_af_noamb_rmN$A1_AF_1000G>0.01) & (b_mega_nodup_af_noamb_rmN$A1_AF_1000G<0.99),]
b_mega_nodup_af_noamb_rmN_comm <- b_mega_nodup_af_noamb_rmN_comm[, c(1:3,6,9:10,12,14,18:20,8)]
colnames(b_mega_nodup_af_noamb_rmN_comm) <- c("rsID", "CHR", "POS_b37", "N", "BETA", "SE", "P", "P_GC", "A1", "A2", "A1_FREQ_1000G", "POOLED_A1_FREQ_glgc")
num_rec[8] <- dim(b_mega_nodup_af_noamb_rmN_comm)[1]; rm("b_mega_nodup_af_noamb_rmN")
write_tsv(b_mega_nodup_af_noamb_rmN_comm, paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",race,"/",trait,".txt"))

names(num_rec) <- c("original","mega","non_missing_p+SNP","remove_duplicated_rsid","in1000G","remove_anbiguous_strand","remove_small_sample_size","common1%")
saveRDS(num_rec, paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",race,"/",trait,"_NofSNP.rds"))


