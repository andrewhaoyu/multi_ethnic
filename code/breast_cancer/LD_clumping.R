#l represent trait
#i1=1 means all SNPs; i1 = 2 means mega snps
#
args = commandArgs(trailingOnly = T)
l = as.numeric(args[[1]])
i1 = as.numeric(args[[2]])
library(data.table)
library(tidyverse)
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')

system(paste0("cp /data/zhangh24/multi_ethnic/data/GBHS_plink/all_chr.bed ",temp.dir,"all_chr.bed"))
system(paste0("cp /data/zhangh24/multi_ethnic/data/GBHS_plink/all_chr.bim ",temp.dir,"all_chr.bim"))
system(paste0("cp /data/zhangh24/multi_ethnic/data/GBHS_plink/all_chr.fam ",temp.dir,"all_chr.fam"))
trait = c("overall","erpos","erneg")
setwd("/data/zhangh24/multi_ethnic/data/")
if(i1 ==1){
  load(paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS.rdata"))
}else{
  load(paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS_mega.rdata"))
}
sum.data.assoc = sum.data %>% 
  select(CHR,ID,POS,Eff_allele,BETA,P) %>% 
  rename(SNP=ID,
         A1 = Eff_allele)
write.table(sum.data.assoc,file = paste0(temp.dir,"sum.data.assoc"),row.names = F,col.names = T,quote=F)
# sum.data.filter = sum.data %>% 
#   filter(Rsq_ave>=0.3&MAF>=0.01)
if(i1==1){
  pthr = 1
  r2thr = 0.1
  kbpthr = 500
  res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",temp.dir,"all_chr --clump ",temp.dir,"sum.data.assoc --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",temp.dir,"LD_clump_",trait[l]))
  system(paste0("mv ",temp.dir,"LD_clump_",trait[l],".clumped /data/zhangh24/multi_ethnic/result/breast_cancer/result/"))
}else{
  pthr = 1
  r2thr = 0.1
  kbpthr = 500
  res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",temp.dir,"all_chr --clump ",temp.dir,"sum.data.assoc --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",temp.dir,"LD_clump_",trait[l],"_mega"))
  system(paste0("mv ",temp.dir,"LD_clump_",trait[l],"_mega.clumped /data/zhangh24/multi_ethnic/result/breast_cancer/result/"))
}
