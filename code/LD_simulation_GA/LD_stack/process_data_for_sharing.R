
#load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/snp.infor.rdata")
# dup.id.list = list()
# library(data.table)
# for(i in 1:5){
#   data = fread(paste0(cur.dir,eth[i],"/duplicated.id"),header=F)
#   dup.id.list[[i]] = data
#   
# }
# 
# dup.id = unique(rbindlist(dup.id.list))
# colnames(dup.id) = "SNP"
#write.table(dup.id,file = paste0(cur.dir,"dup.id"),col.names = T,row.names = F,quote=F)
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
# args = commandArgs(trailingOnly = T)
# i = as.numeric(args[[1]])
# j = as.numeric(args[[2]])
# eth = c("EUR","AFR","AMR","EAS","SAS")
# 
# 
# library(data.table)
# library(dplyr)
# sid<-Sys.getenv('SLURM_JOB_ID')
# dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
# temp.dir = paste0('/lscratch/',sid,'/test/')
# cur.dir = "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
# out.dir <-  paste0(cur.dir,eth[i],"/",eth[i],"_mega/")
# 
# system(paste0("cp ",cur.dir,eth[i],"/chr",j,".mega.bed ",temp.dir,eth[i],"chr",j,".mega.bed"))
# system(paste0("cp ",cur.dir,eth[i],"/chr",j,".mega.bim ",temp.dir,eth[i],"chr",j,".mega.bim"))
# system(paste0("cp ",cur.dir,eth[i],"/chr",j,".mega.fam ",temp.dir,eth[i],"chr",j,".mega.fam"))
# 
# system(paste0("cp ",cur.dir,eth[i],"/chr",j,"train.mega.bed ",temp.dir,eth[i],"chrtrain",j,".mega.bed"))
# system(paste0("cp ",cur.dir,eth[i],"/chr",j,"train.mega.bim ",temp.dir,eth[i],"chrtrain",j,".mega.bim"))
# system(paste0("cp ",cur.dir,eth[i],"/chr",j,"train.mega.fam ",temp.dir,eth[i],"chrtrain",j,".mega.fam"))
# 
# all_my_files <- rep("c",1)
#   all_my_files[1] = paste0(temp.dir,eth[i],"chrtrain",j,".mega")
# write.table(all_my_files,file = paste0(temp.dir,"all_my_files.txt"),row.names = F,col.names = F,quote=F)
# 
# res = system(paste0("/data/zhangh24/software/plink2 ",
#                     "--bfile ",temp.dir,eth[i],"chr",j,".mega ",
#                     "--merge-list ",temp.dir,"all_my_files.txt ",
#                     "--out ",temp.dir,"chr",j," ",
#                     "--exclude ",cur.dir,"dup.id ",
#                     "--make-bed "))
# 
# system(paste0("cd ",temp.dir,"; zip -r chr",j,".bed.zip chr",j,".bed"))
# system(paste0("cd ",temp.dir,"; zip -r chr",j,".bim.zip chr",j,".bim"))
# system(paste0("cd ",temp.dir,"; zip -r chr",j,".fam.zip chr",j,".fam"))
# system(paste0("cp ",temp.dir,"chr",j,".bed.zip ",out.dir,"chr",j,".bed.zip"))
# system(paste0("cp ",temp.dir,"chr",j,".fam.zip ",out.dir,"chr",j,".fam.zip"))
# system(paste0("cp ",temp.dir,"chr",j,".bim.zip ",out.dir,"chr",j,".bim.zip"))







args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
eth = c("EUR","AFR","AMR","EAS","SAS")


library(data.table)
library(dplyr)
cur.dir = "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir <-  paste0(cur.dir,eth[i],"/",eth[i],"_mega/")

system(paste0("cd " ,cur.dir,"; ",
"zip -r ",eth[i],"_mega.zip ",eth[i],"_mega/"))

for(j in 1:22){
  system(paste0("cd ",cur.dir,eth[i],"; ",
                "rm chr",j,"train.mega.* ; ",
                "rm chr",j,".mega.* "))
}



