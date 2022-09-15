#goal: create the bim file for prs-csx five ancestries
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"
library(data.table)
library(dplyr)
snp.list = list()
l = 1
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")
data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"


#the SNP ID for different traits are consistent
#only need to create a single bim for all data analyses
for(i in 1:5){
  sum = as.data.frame(fread(paste0(data.dir,eth[i],"/sumdat/",trait[3],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
  bim1 = sum  %>% 
    mutate(V3=0) %>% 
    select(CHR,rsid,V3,BP,A1,A2)
  
  sum = as.data.frame(fread(paste0(data.dir,eth[i],"/sumdat/",trait[4],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
  bim2 = sum  %>% 
    mutate(V3=0) %>% 
    select(CHR,rsid,V3,BP,A1,A2)
  print(all.equal(bim1,bim2))
  
}


for(i in 1:5){
  sum = as.data.frame(fread(paste0(data.dir,eth[i],"/sumdat/",trait[2],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
  
  snp.list[[i]] = sum  %>% 
    mutate(V3=0) %>% 
    select(CHR,rsid,V3,BP,A1,A2)
  
}


snp.file = rbindlist(snp.list) %>% 
  distinct(rsid,.keep_all=TRUE)
write.table(snp.file,file = paste0(data.dir,"all_eth.bim"),row.names = F,col.names = F,quote=F)
