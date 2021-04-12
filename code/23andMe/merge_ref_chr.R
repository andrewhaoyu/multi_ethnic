#goal: merge ARIC different chrs into one dataset
eth <- c("EUR","AFR","AMR","EAS","SAS")


#i = 1
eth <- c("EUR","AFR","AMR","EAS","SAS")
for(i in 1:5){
  data.dir = paste0("/data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"/")
  temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
  
  total = 21
  temp = 1
  filename = rep("c",total)
  for(j in 2:22){
    filename[temp] = paste0(data.dir,"chr",j)
    temp = temp+1
  }
  out.dir = data.dir
  write.table(filename,file = paste0(out.dir,"merge_list.txt"),row.names = F,col.names = F,quote=F)
  
  system(paste0("/data/zhangh24/software/plink2 --bfile ",data.dir,"chr",1," --merge-list ",data.dir,"merge_list.txt --make-bed --out ",data.dir,"all_chr"))
  
}

#transform into rds using bigsnpr
# library(bigsnpr)
# snp_readBed(paste0(temp.dir,"all_chr.bed"))


data <- fread("/gpfs/gsfs11/users/zhangh24/KGref_MEGA/GRCh37/EUR/all_chr.bim")

idx <- which(data$V2%in%c("rs201335322","rs10539207"))
data[idx,]
