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
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
#load EUR clumping summary data
summary.eur <- as.data.frame(fread(paste0(out.dir.sum,eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))  
colnames(summary.eur)[9] = "peur"
summary.eur = summary.eur %>% filter(CHR==22)

#for(i in 2:4){
#summary <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/summary.out"),header=T))
#load target snps summary stat
summary <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
summary.tar = summary %>% filter(CHR==22)
#load 1kg snp list
setwd("/data/zhangh24/multi_ethnic/")
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
snp.infor = snp.infor.match %>% 
  rename(SNP=id) %>% 
  select(SNP,rs_id)
summary.eur = left_join(summary.eur,snp.infor,by="SNP")
summary.tar = left_join(summary.tar,snp.infor,by="SNP")
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

#prepare the data for clumping reference
snp.list = unique(c(summary.eur$SNP,summary.tar$SNP))
write.table(snp.list,file = paste0(out.dir.sum,"time_mem/extract_snp_list.txt"),row.names = F,col.names = F,quote=F)
system(paste0("cp ",cur.dir,eth[1],"/clump_ref_all_chr.bed ",temp.dir,eth[1],"clump_ref_all_chr.bed"))
res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",cur.dir,eth[1],"/clump_ref_all_chr --extract ",out.dir.sum,"time_mem/extract_snp_list.txt --out ", out.dir.sum,"time_mem/EUR_ref_chr22 --make-bed"))
res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",cur.dir,eth[2],"/clump_ref_all_chr --extract ",out.dir.sum,"time_mem/extract_snp_list.txt --out ", out.dir.sum,"time_mem/AFR_ref_chr22 --make-bed"))

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
