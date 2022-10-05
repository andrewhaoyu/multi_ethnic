#test the running time and memory for TDLD-SLEB and PRS-CSx on CHR 22
i = 2
l = 3
m = 1
i_rep = 1
i1 = 1
library(data.table)
library(dplyr)
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/LD/'),showWarnings = FALSE)
temp.dir.LD <- paste0('/lscratch/',sid,'/test/LD/')
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/time_mem/"
#load EUR clumping summary data
summary.eur <- as.data.frame(fread(paste0(out.dir.sum,eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))  
colnames(summary.eur)[9] = "peur"
summary.eur = summary.eur %>% filter(CHR==22)

#for(i in 2:4){
#summary <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out"),header=T))
#load target snps summary stat
summary <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
summary.tar = summary %>% filter(CHR==22)

summary.amr <- as.data.frame(fread(paste0(out.dir.sum,eth[3],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
summary.amr = summary.amr %>% filter(CHR==22)
summary.eas <- as.data.frame(fread(paste0(out.dir.sum,eth[4],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
summary.eas = summary.eas %>% filter(CHR==22)
summary.sas <- as.data.frame(fread(paste0(out.dir.sum,eth[5],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
summary.sas = summary.sas %>% filter(CHR==22)







#load 1kg snp list
setwd("/data/zhangh24/multi_ethnic/")
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
snp.infor = snp.infor.match %>% 
  rename(SNP=id) %>% 
  select(SNP,rs_id)
summary.eur = left_join(summary.eur,snp.infor,by="SNP")
summary.tar = left_join(summary.tar,snp.infor,by="SNP")
summary.amr = left_join(summary.amr,snp.infor,by="SNP")
summary.eas = left_join(summary.eas,snp.infor,by="SNP")
summary.sas = left_join(summary.sas,snp.infor,by="SNP")
mega.list <- as.data.frame(fread(paste0(cur.dir,"mega-hm3-rsid.txt"),header  =F))
#mega.infor <- as.data.frame(fread(paste0(cur.dir,"snpBatch_ILLUMINA_1062317")))
#colnames(mega.infor)[5] <- "rsid"
colnames(mega.list) = "rs_id"

#subset SNPs to mega SNPs
summary.eur.match = inner_join(summary.eur,mega.list,by="rs_id")
summary.tar.match = inner_join(summary.tar,mega.list,by="rs_id")
summary.tar.match = summary.tar.match %>% 
  select(CHR,SNP,BP,A1,TEST,NMISS,BETA,STAT,P,rs_id)
summary.eur = summary.eur.match
summary.tar = summary.tar.match
save(summary.eur,file = paste0(out.dir.sum,"time_mem/eur_sumdata.rdata"))
save(summary.tar,file = paste0(out.dir.sum,"time_mem/tar_sumdata.rdata"))
save(summary.amr,file = paste0(out.dir.sum,"time_mem/amr_sumdata.rdata"))
save(summary.eas,file = paste0(out.dir.sum,"time_mem/eas_sumdata.rdata"))
save(summary.sas,file = paste0(out.dir.sum,"time_mem/sas_sumdata.rdata"))
#prepare SNPs for other ethnic groups
summary.eur.select = summary.eur %>% 
  rename(beta_eur = BETA) %>% 
  mutate(sd_eur=beta_eur/STAT) %>% 
  select(SNP,A1,beta_eur,sd_eur,peur) %>% 
  rename(A1.EUR = A1)
summary.com <- left_join(summary.tar,summary.eur.select,by="SNP")
i_can <- setdiff(c(2:5),i)
beta.mat.list = list()
sd.mat.list = list()
temp = 1



for(i_sub in i_can){
  sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i_sub],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))
  colnames(sum.data)[2] <- "SNP"
  sum.data = sum.data %>% 
    mutate(sd.sub = BETA/STAT) %>% 
    rename(A1.sub = A1,
           BETA.sub = BETA) %>% 
    select(SNP,A1.sub,BETA.sub,sd.sub)
  summary.com.temp <- left_join(summary.com,sum.data,by="SNP")
  #match allele
  summary.com.temp = summary.com.temp %>% 
    mutate(BETA.sub.new = ifelse(A1==A1.sub,BETA.sub,-BETA.sub)) 
  beta.mat.list[[temp]] =  summary.com.temp %>% select(BETA.sub.new)
  sd.mat.list[[temp]] = summary.com.temp %>% select(sd.sub)
  temp = temp+1
}
beta.mat = bind_cols(beta.mat.list)
sd.mat = bind_cols(sd.mat.list)
save(beta.mat,file=paste0(out.dir, "other_ethnic_beta_mat"))
save(sd.mat,file = paste0(out.dir,"other_ethnic_sd_mat"))
#prepare the data for clumping reference
snp.list = unique(c(summary.eur$SNP,summary.tar$SNP))
write.table(snp.list,file = paste0(out.dir.sum,"time_mem/extract_snp_list.txt"),row.names = F,col.names = F,quote=F)
system(paste0("cp ",cur.dir,eth[1],"/clump_ref_all_chr.bed ",temp.dir,eth[1],"clump_ref_all_chr.bed"))
res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",cur.dir,eth[1],"/clump_ref_all_chr --extract ",out.dir.sum,"time_mem/extract_snp_list.txt --out ", out.dir.sum,"time_mem/EUR_ref_chr22 --make-bed"))
res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",cur.dir,eth[2],"/clump_ref_all_chr --extract ",out.dir.sum,"time_mem/extract_snp_list.txt --out ", out.dir.sum,"time_mem/AFR_ref_chr22 --make-bed"))


# system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile /data/zhangh24/KG.plink/EUR/chr_all --extract ",out.dir.sum,"time_mem/extract_snp_list.txt --out /data/zhangh24/KG.plink/EUR/EUR_ref_1kg_chr22 --make-bed"))
# system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile /data/zhangh24/KG.plink/AFR/chr_all --extract ",out.dir.sum,"time_mem/extract_snp_list.txt --out /data/zhangh24/KG.plink/AFR/AFR_ref_1kg_chr22 --make-bed"))




#prepare the data for prs calculation
system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",cur.dir,eth[2],"/all_chr_test.mega --extract ",out.dir.sum,"time_mem/extract_snp_list.txt --out ", out.dir.sum,"time_mem/AFR_test_mega_chr22 --make-bed"))
              
#prepare phenotype data for r2 calculation
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
y <- as.data.frame(fread(paste0(cur.dir,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
n.rep = 10
y <- y[,2+(1:n.rep),drop=F]
n <- nrow(y)
y_test_mat <- y[(100000+1):nrow(y),,drop=F]
y_test = y_test_mat[1:10000,i_rep]
y_vad = y_test_mat[10001:20000,i_rep]
save(y_test,file = paste0(out.dir,"y_test.rdata"))
save(y_vad,file = paste0(out.dir,"y_vad.rdata"))

#prepare data for prscsx
setwd("/data/zhangh24/multi_ethnic/")
load(paste0(out.dir.sum,"time_mem/eur_sumdata.rdata"))
prs_cs_ref = read.table("/gpfs/gsfs11/users/zhangh24/software/PRScsx/1KGLD/snpinfo_mult_1kg_hm3",header=T)
prs_cs_ref = prs_cs_ref %>% select(SNP)
summary.eur = inner_join(summary.eur,prs_cs_ref,by = c("rs_id"="SNP"))
library(tidyverse)
summary.eur = summary.eur %>% 
  separate(SNP,into=c("non_im","non_im2","first_allele","second_allele"),
           remove=F) %>% 
  mutate(A2 =ifelse(A1==second_allele,first_allele,second_allele)) %>% 
  rename(P = peur) %>% 
  select(rs_id,A1,A2,BETA,P) %>% 
  rename(SNP=rs_id)
write.table(summary.eur,file = paste0(out.dir,"EUR_sumstats.txt"),row.names = F,col.names = T,quote=F)
#merge with prs-csx reference data
#prepare target summary stat for a specific chr
load(paste0(out.dir.sum,"time_mem/tar_sumdata.rdata"))
#merge with prs-csx reference data
summary.tar = inner_join(summary.tar,prs_cs_ref,by = c("rs_id"="SNP"))
summary.tar.sub = summary.tar %>% 
  separate(SNP,into=c("non_im","non_im2","first_allele","second_allele"),
           remove=F) %>% 
  mutate(A2 =ifelse(A1==second_allele,first_allele,second_allele)) %>% 
  select(rs_id,A1,A2,BETA,P) %>% 
  rename(SNP=rs_id)
write.table(summary.tar.sub,file = paste0(out.dir,eth[i],"_sumstats.txt"),row.names = F,col.names = T,quote=F)
#prepare amr summary stat for a specific chr
load(paste0(out.dir.sum,"time_mem/amr_sumdata.rdata"))
#merge with prs-csx reference data
summary.amr = inner_join(summary.amr,prs_cs_ref,by = c("rs_id"="SNP"))
summary.amr.sub = summary.amr %>% 
  separate(SNP,into=c("non_im","non_im2","first_allele","second_allele"),
           remove=F) %>% 
  mutate(A2 =ifelse(A1==second_allele,first_allele,second_allele)) %>% 
  select(rs_id,A1,A2,BETA,P) %>% 
  rename(SNP=rs_id)
write.table(summary.amr.sub,file = paste0(out.dir,"AMR_sumstats.txt"),row.names = F,col.names = T,quote=F)

#prepare EAS summary stat for a specific chr
load(paste0(out.dir.sum,"time_mem/eas_sumdata.rdata"))
#merge with prs-csx reference data
summary.eas = inner_join(summary.eas,prs_cs_ref,by = c("rs_id"="SNP"))
summary.eas.sub = summary.eas %>% 
  separate(SNP,into=c("non_im","non_im2","first_allele","second_allele"),
           remove=F) %>% 
  mutate(A2 =ifelse(A1==second_allele,first_allele,second_allele)) %>% 
  select(rs_id,A1,A2,BETA,P) %>% 
  rename(SNP=rs_id)
write.table(summary.eas.sub,file = paste0(out.dir,"EAS_sumstats.txt"),row.names = F,col.names = T,quote=F)

#prepare SAS summary stat for a specific chr
load(paste0(out.dir.sum,"time_mem/sas_sumdata.rdata"))
#merge with prs-csx reference data
summary.sas = inner_join(summary.sas,prs_cs_ref,by = c("rs_id"="SNP"))
summary.sas.sub = summary.sas %>% 
  separate(SNP,into=c("non_im","non_im2","first_allele","second_allele"),
           remove=F) %>% 
  mutate(A2 =ifelse(A1==second_allele,first_allele,second_allele)) %>% 
  select(rs_id,A1,A2,BETA,P) %>% 
  rename(SNP=rs_id)
write.table(summary.eas.sub,file = paste0(out.dir,"SAS_sumstats.txt"),row.names = F,col.names = T,quote=F)



bim = read.table(paste0(out.dir,"AFR_ref_chr22.bim"))
bim.update =inner_join(summary.tar,bim,by = c("SNP"="V2")) %>% 
  select(V1,rs_id,V3,V4,V5,V6)
write.table(bim.update,file = paste0(out.dir,"AFR_ref_chr22_rs_id.bim"),row.names = F,col.names = F,quote=F)
#phi = c(1E-6,1E-4)

fam = fread(paste0(out.dir.sum,"/time_mem/AFR_test_mega_chr22.fam"),header=F)
fam_test = fam[1:10000,]
fam_vad = fam[1:10000,]
write.table(fam_test, file = paste0(out.dir.sum,"/time_mem/test_id.fam"),row.names = F,col.names = F,quote=F)
write.table(fam_vad, file = paste0(out.dir.sum,"/time_mem/vad_id.fam"),row.names = F,col.names = F,quote=F)
#prepare the data for CTSLEB demonstration
system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",out.dir.sum,"time_mem/AFR_test_mega_chr22 --keep ",out.dir.sum,"time_mem/test_id.fam --out ", out.dir.sum,"time_mem/AFR_test_chr22 --make-bed"))
system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",out.dir.sum,"time_mem/AFR_test_mega_chr22 --keep ",out.dir.sum,"time_mem/vad_id.fam --out ", out.dir.sum,"time_mem/AFR_vad_chr22 --make-bed"))
