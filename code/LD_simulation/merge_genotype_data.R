args <- commandArgs(trailingOnly = T)
i <- as.numeric(args[[1]])
j <- as.numeric(args[[2]])
#number of jobs
n <- 1000
eth <- c("EUR","AFR","AMR","EAS","SAS")
merge_list = rep("c",n-1)
temp = 1
#update the family id and individual id in the fam file
library(data.table)
total <- 0
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
for(k in 1:n){
  fam <- read.table(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i],"/chr",j,"_",k,".tag.fam"),header=F)
  temp <- nrow(fam)
  fam[,1] <- c(total+(1:temp))
  fam[,2] <- c(total+(1:temp))
  write.table(fam,file= paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i],"/chr",j,"_",k,".tag.fam"),
              col.names = F,row.names = F,quote=F)
  total <- total+temp
}
temp = 1
for(k in 2:n){
  merge_list[temp] <- paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i],"/chr",j,"_",k,".tag")
  temp = temp+1
}
write.table(merge_list,file = paste0("/lscratch/",sid,"/test/",eth[i],"_all_files_chr_",j),
            col.names = F,
            row.names = F,
            quote=F)
#merge all the seperate sample into one
#merge the seperate genotype data into one file
system(paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/multi_ethnic/result/LD_simulation_new/EUR/chr2_1.tag --merge-list /lscratch/",sid,"/test/",eth[i],"_all_files_chr_",j, " --make-bed --out /data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i],"/chr",j,".tag"))






