args = commandArgs(trailingOnly = T)
l = as.numeric(args[[1]])
i1 = as.numeric(args[[2]])
library(tidyverse)
library(data.table)
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5,1.0)
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
out.dir = "/data/zhangh24/multi_ethnic/result/breast_cancer/result/"
system(paste0("cp /data/zhangh24/multi_ethnic/data/GBHS_plink/all_chr.bed ",temp.dir,"all_chr.bed"))
system(paste0("cp /data/zhangh24/multi_ethnic/data/GBHS_plink/all_chr.bim ",temp.dir,"all_chr.bim"))
system(paste0("cp /data/zhangh24/multi_ethnic/data/GBHS_plink/all_chr.fam ",temp.dir,"all_chr.fam"))
trait = c("overall","erpos","erneg")
setwd("/data/zhangh24/multi_ethnic/data/")

out.dir <-  "/data/zhangh24/multi_ethnic/result/breast_cancer/result/"
if(i1 ==1){
  load(paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS.rdata"))
  LD <- as.data.frame(fread(paste0(out.dir,"/LD_clump_",trait[l],".clumped")))
}else{
  load(paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS_mega.rdata"))
  LD <- as.data.frame(fread(paste0(out.dir,"/LD_clump_",trait[l],"_mega.clumped")))
}




clump.snp <- LD[,3,drop=F] 
sum.data = sum.data %>% 
  select(CHR,ID,POS,Eff_allele,BETA,P) %>% 
  rename(SNP=ID,
         A1 = Eff_allele)

prs.all <- left_join(clump.snp,sum.data,by="SNP")
#find all the duplicated SNPs and keep the more significant one
dup.id <- prs.all$SNP[duplicated(prs.all$SNP)]
remove.idx = rep(0,length(dup.id))
for(k in 1:length(dup.id)){
  jdx <- which(prs.all$SNP==dup.id[k])
  which.max(prs.all$P[jdx])
  remove.idx[k] = jdx[which.max(prs.all$P[jdx])]
}
prs.all = prs.all[-remove.idx,]
n_pthres <- length(pthres)
q_range = data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres))

prs.file <- prs.all %>% 
  select(SNP,A1,BETA)
write.table(prs.file,file = paste0(temp.dir,"prs_file"),col.names = T,row.names = F,quote=F)
p.value.file <- prs.all %>% 
  select(SNP,P)
write.table(p.value.file,file = paste0(temp.dir,"p_value_file"),col.names = T,row.names = F,quote=F)
if(i1==1){
  res <- system(paste0("/data/zhangh24/software/plink2 --q-score-range ",temp.dir,"q_range_file ",temp.dir,"p_value_file header --threads 2 --score ",temp.dir,"prs_file header no-sum no-mean-imputation --bfile ",temp.dir,"all_chr --out ",temp.dir,"prs_",trait[l]))
  res = system(paste0("mv ",temp.dir,"*.profile ",out.dir,))
}else{
  res <- system(paste0("/data/zhangh24/software/plink2 --q-score-range ",temp.dir,"q_range_file ",temp.dir,"p_value_file header --threads 2 --score ",temp.dir,"prs_file header no-sum no-mean-imputation --bfile ",temp.dir,"all_chr --out ",temp.dir,"prs_",trait[l],"_mega"))
}

