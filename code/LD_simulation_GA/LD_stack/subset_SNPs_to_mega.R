#goal subset the SNPs data to mega list
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])
i1 = as.numeric(args[[5]])
library(data.table)
library(dplyr)

sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
system(paste0("cp ",cur.dir,eth[i],"/chr",j,".tag.bed ",temp.dir,eth[i],"chr",j,".tag.bed"))
system(paste0("cp ",cur.dir,eth[i],"/chr",j,".tag.bim ",temp.dir,eth[i],"chr",j,".tag.bim"))
system(paste0("cp ",cur.dir,eth[i],"/chr",j,".tag.fam ",temp.dir,eth[i],"chr",j,".tag.fam"))

MAF.cutoff = 0.005
if(i==1){
  MAF.cutoff = 0.01
}
library(rlang)
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
snp.infor = snp.infor.match %>% 
  filter(!!sym(eth[i])>=MAF.cutoff&
           !!sym(eth[i])<=(1-MAF.cutoff)&
           CHR==j) %>% 
  rename(SNP=id) %>% 
  select(SNP,rs_id)

mega.list <- as.data.frame(fread(paste0(cur.dir,"mega-hm3-rsid.txt"),header  =F))
#mega.infor <- as.data.frame(fread(paste0(cur.dir,"snpBatch_ILLUMINA_1062317")))
#colnames(mega.infor)[5] <- "rsid"
colnames(mega.list) = "rs_id"

snp.infor.subset = inner_join(snp.infor,mega.list,by="rs_id") %>% 
  select(SNP)
write.table(snp.infor.subset,file = paste0(temp.dir,"extract_snp_list.txt"),row.names = F,col.names = F,quote=F)

res = system(paste0("/data/zhangh24/software/plink2 --bfile ",temp.dir,eth[i],"chr",j,".tag --extract ",temp.dir,"extract_snp_list.txt --out ",temp.dir,"chr",j,".mega --make-bed"))

if(res==2){
  stop()
}

system(paste0("mv ",temp.dir,"/chr",j,".mega.bed ",cur.dir,eth[i],"/chr",j,".mega.bed"))
system(paste0("mv ",temp.dir,"/chr",j,".mega.bim ",cur.dir,eth[i],"/chr",j,".mega.bim"))
system(paste0("mv ",temp.dir,"/chr",j,".mega.fam ",cur.dir,eth[i],"/chr",j,".mega.fam"))
