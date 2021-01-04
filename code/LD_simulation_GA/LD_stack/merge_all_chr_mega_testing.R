#goal subset the SNPs data to mega list
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
library(data.table)
library(dplyr)

sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test/'))
temp.dir = paste0('/lscratch/',sid,'/test/')
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
for(j in 1:22){
  system(paste0("cp ",cur.dir,eth[i],"/chr",j,".mega.bed ",temp.dir,eth[i],"chr",j,".mega.bed"))
  system(paste0("cp ",cur.dir,eth[i],"/chr",j,".mega.bim ",temp.dir,eth[i],"chr",j,".mega.bim"))
  system(paste0("cp ",cur.dir,eth[i],"/chr",j,".mega.fam ",temp.dir,eth[i],"chr",j,".mega.fam"))
}

# MAF.cutoff = 0.005
# if(i==1){
#   MAF.cutoff = 0.01
# }
library(rlang)
setwd("/data/zhangh24/multi_ethnic/")
# load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
# snp.infor = snp.infor.match %>% 
#   filter(!!sym(eth[i])>=MAF.cutoff&
#            !!sym(eth[i])<=(1-MAF.cutoff)&
#            CHR==j) %>% 
#   rename(SNP=id) %>% 
#   select(SNP,rs_id)
# 
# mega.list <- as.data.frame(fread(paste0(cur.dir,"mega-hm3-rsid.txt"),header  =F))
# #mega.infor <- as.data.frame(fread(paste0(cur.dir,"snpBatch_ILLUMINA_1062317")))
# #colnames(mega.infor)[5] <- "rsid"
# colnames(mega.list) = "rs_id"
# 
# snp.infor.subset = inner_join(snp.infor,mega.list,by="rs_id") %>% 
#   select(SNP)
# write.table(snp.infor.subset,file = paste0(temp.dir,"extract_snp_list.txt"),row.names = F,col.names = F,quote=F)
all_my_files <- rep("c",21)
temp = 1
for(j in 2:22){
  all_my_files[temp] = paste0(temp.dir,eth[i],"chr",j,".mega")
  temp = temp+1
  
}
write.table(all_my_files,file = paste0(temp.dir,"all_my_files.txt"),row.names = F,col.names = F,quote=F)
res = system(paste0("/data/zhangh24/software/plink2 --bfile ",temp.dir,eth[i],"chr",1,".mega --merge-list ",temp.dir,"all_my_files.txt --out ",temp.dir,"all_chr_test.mega --make-bed"))

if(res==2){
  stop()
}

system(paste0("mv ",temp.dir,"/all_chr_test.mega.bed ",cur.dir,eth[i],"/all_chr_test.mega.bed"))
system(paste0("mv ",temp.dir,"/all_chr_test.mega.bim ",cur.dir,eth[i],"/all_chr_test.mega.bim"))
system(paste0("mv ",temp.dir,"/all_chr_test.mega.fam ",cur.dir,eth[i],"/all_chr_test.mega.fam"))
system(paste0("rm -rf /lscratch/",sid,"/test"))
system(paste0("ls /lscratch/",sid))
