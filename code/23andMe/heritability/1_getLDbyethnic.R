#calculate the LD-score for each ethnic using 1000 genome data
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
j = as.numeric(args[[2]])
eth = c("EUR","AFR","AMR","EAS","SAS")
#create temp file folder
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/',eth[i],"/")
system(paste0("mkdir ",temp.dir))
kg.dir = "/data/zhangh24/KGref_MEGA/GRCh37/"

system(paste0("cp ",kg.dir,eth[i],"/chr",j,".bed ",temp.dir,"chr",j,".bed"))
system(paste0("cp ",kg.dir,eth[i],"/chr",j,".bim ",temp.dir,"chr",j,".bim"))
system(paste0("cp ",kg.dir,eth[i],"/chr",j,".fam ",temp.dir,"chr",j,".fam"))

#remove HLC region (25002566-34994414)
library(dplyr)
if(j==6){
  bim = read.table(paste0(temp.dir,"chr",j,".bim"))
  head(bim)
  bim.hla = bim %>% 
    filter(V4>=25002566&
             V4<=34994414)
  dim(bim.hla)
  write.table(bim.hla$V2,file = paste0(temp.dir,"remove_region"),row.names = F,col.names = F,quote=F)
  system(paste0("/data/zhangh24/software/plink2 --bfile ",temp.dir,"chr",j," --exclude ",temp.dir,"remove_region --make-bed --out ",temp.dir,"chr_clean",j))
  system(paste0("source activate ldsc; ", 
                "cd /data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"; ",
                "python /data/zhangh24/ldsc/ldsc.py ",
                "--bfile ",temp.dir,"chr_clean",j," ",
                "--l2 --ld-wind-kb 1000 ",
                "--out /data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"/",eth[i],"_ldsc/",j))
  
}else{
  system(paste0("source activate ldsc; ", 
                "cd /data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"; ",
                "python /data/zhangh24/ldsc/ldsc.py ",
                "--bfile ",temp.dir,"chr",j," ",
                "--l2 --ld-wind-kb 1000 ",
                "--out /data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"/",eth[i],"_ldsc/",j))
  
}

