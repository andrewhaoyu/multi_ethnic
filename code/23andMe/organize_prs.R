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
total.prs = length(pthres)
#load snp_info file to get corresponding im.id
load(paste0(data.dir,"snpinfo/snpinfo_mega.RData"))
snpinfo_mega_filter = snpinfo_mega %>% 
  filter(!is.na(im.data.id)) %>% 
  select(im.data.id,assay.name)
#genearte index file to match columns in organize prs and parameters
total = length(eth)*length(trait)*length(pthres)
eth_vec = rep("c",total)
trait_vec = rep("c",total)
pvalue_vec = rep(0,total)
col_vec = rep(0,total)
temp = 1
for(i in 1:length(eth)){
  for(l in 1:length(trait)){
    out.dir.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/PT/",eth[i],"/",trait[l],"/")
    files = dir(out.dir.prs,method)
    files.temp = gsub("_PT_pvalue","",files)
    files.split = strsplit(files.temp,"_")
    for(k in 1:length(files)){
      col_vec[temp] = as.numeric(files.split[[k]][1])
      pvalue_vec[temp] = as.numeric(files.split[[k]][2])
      eth_vec[temp] = eth[i]
      trait_vec[temp] = trait[l]
      temp =temp+1
    }
  }
  }
total = temp-1
eth_vec = eth_vec[1:total]
trait_vec = trait_vec[1:total]
pvalue_vec = pvalue_vec[1:total]
col_vec = col_vec[1:total]

index.match = data.frame(eth_vec,trait_vec,pvalue_vec,col_vec)
save(index.match,file = "/data/zhangh24/multi_ethnic/result/cleaned/prs/PT/index.match.rdata")

for(i in 1:length(eth)){
  for(l in 1:length(trait)){
    print(c(i,l))
    out.dir.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/PT/",eth[i],"/",trait[l],"/")
    temp = 1
    rs.id.list = list()
    #run through all the prs to get unique number of snps
    files = dir(out.dir.prs,method)
    for(k in 1:length(files)){
      prs.file <- fread(paste0(out.dir.prs,files[k]),header=T)
      rs.id.list[[temp]] = prs.file[,1,drop=F]
      temp = temp+1
    }
    rs.id = unique(rbindlist(rs.id.list))
    
    #prs snp list
    #order snps by im.data.id
    prs.snp = left_join(rs.id,snpinfo_mega_filter,by=c("SNP"="assay.name")) %>% 
      arrange(im.data.id)
    temp = 1
    files = dir(out.dir.prs,method)
    for(k in 1:length(files)){
      prs.file <- fread(paste0(out.dir.prs,files[k]),header=T)
      prs.snp = left_join(prs.snp,prs.file,by="SNP") %>% 
        mutate(BETA = ifelse(is.na(BETA),0,BETA))
      colnames(prs.snp)[temp+2] = paste0("BETA",temp)
      temp = temp+1
    }
    prs.snp = prs.snp[,-1]
    out.dir.organize.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/organize_prs/PT/",eth[i],"/",trait[l],"/")
    write.table(prs.snp,file = paste0(out.dir.organize.prs,"prs.file"),row.names = F,col.names = F,quote=F)
    
  }
}