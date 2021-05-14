args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)
method = "PT"
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"
#load best eur results
out.dir = "/data/zhangh24/multi_ethnic/result/cleaned/result_summary/"
eth_group = c("european","african_american",
              "latino","east_asian","south_asian")
method = "PT"
result <- read.csv(paste0(out.dir,eth_group[1],"_",trait[l],"_",method))
k <- which.max(result)
out.dir.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/PT/",eth[1],"/",trait[l],"/")
prs.file.name = dir(out.dir.prs,pattern = paste0("_",method,"_pvalue_",k),full.names = T)
#best eur prs
best.prs = read.table(prs.file.name,header=T)
best.eur.prs = best.prs %>% 
  rename(BETA.EUR = BETA)
#load EUR ethnic group coefficients
sum.data = as.data.frame(fread(paste0(data.dir,eth[1],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
best.eur.prs = left_join(best.eur.prs,sum.data,by=c("SNP"="rsid")) %>% 
  select(SNP,FREQ_A1,BETA,SD) %>% 
  rename(BETA.EUR = BETA,
         SD.EUR = SD,
         FREQ.EUR = FREQ_A1)

#load target ethnic group coefficients
sum.data = as.data.frame(fread(paste0(data.dir,eth[i],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
#find the subset of best eur prs in target ethnic group
best.eur.prs.com = inner_join(best.eur.prs,sum.data,by=c("SNP"="rsid")) %>% 
  select(SNP,BETA.EUR,SD.EUR,FREQ.EUR,BETA,SD,FREQ_A1) %>% 
  rename(BETA.TAR = BETA,
         SD.TAR = SD,
         FREQ.TAR = FREQ_A1)
source("/data/zhangh24/multi_ethnic/code/stratch/EB_function.R")

beta_tar <- best.eur.prs.com$BETA.TAR
sd_tar <- best.eur.prs.com$SD.TAR
beta_eur <- best.eur.prs.com$BETA.EUR
sd_eur <- best.eur.prs.com$SD.EUR

EBprior = EstimatePrior(beta_tar,sd_tar,
                        beta_eur,sd_eur)

post_beta_mat = EBpost(beta_tar,sd_tar,beta_eur,sd_eur,EBprior)
post_beta_tar = post_beta_mat[,1,drop=F]

best.eur.prs.com$BETA.EB = post_beta_tar


#best target ethnic group prs results
method = "PT"
result <- read.csv(paste0(out.dir,eth_group[i],"_",trait[l],"_",method))
k <- which.max(result)
out.dir.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/PT/",eth[i],"/",trait[l],"/")
prs.file.name = dir(out.dir.prs,pattern = paste0("_",method,"_pvalue_",k),full.names = T)
#best target prs
best.prs = read.table(prs.file.name,header=T)


#combime best eur prs and best target ethnic group prs

best.prs.all = full_join(best.eur.prs.com,
                         best.prs,by="SNP")
#fill in 0 for the rows without any cofficients
best.prs.select = best.prs.all %>% 
  select(SNP,BETA.EUR,BETA.TAR,BETA.EB,BETA)
best.prs.select[is.na(best.prs.select)] = 0
#match snp with imputed id
load(paste0(data.dir,"snpinfo/snpinfo_mega.RData"))
snpinfo_mega_filter = snpinfo_mega %>% 
  filter(!is.na(im.data.id)) %>% 
  select(im.data.id,assay.name)
prs.snp = left_join(best.prs.select,snpinfo_mega_filter,by=c("SNP"="assay.name")) %>% 
  arrange(im.data.id) %>% 
  select(im.data.id,BETA.EUR,BETA.TAR,BETA.EB,BETA)

out.dir.organize.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/organize_prs/BESTEUR/",eth[i],"/",trait[l],"/")
write.table(prs.snp,file = paste0(out.dir.organize.prs,"prs.file"),row.names = F,col.names = F,quote=F)