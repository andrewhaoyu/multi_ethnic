#combine the by chromosome data into one file
args <- commandArgs(trailingOnly = T)
i <- as.numeric(args[[1]])
n <- 22
eth <- c("EUR","AFR","AMR","EAS","SAS")
merge_list = rep("c",21)
temp = 1
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
for(j in 2:n){
  merge_list[temp] <- paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i],"/chr",j,".tag")
  temp = temp+1
}
write.table(merge_list,file = paste0("/lscratch/",sid,"/test/",eth[i],"_all_files"),
            col.names = F,
            row.names = F,
            quote=F)
#merge all the seperate sample into one
#merge the seperate genotype data into one file
system(paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i],"/chr",1,".tag --merge-list /lscratch/",sid,"/test/",eth[i],"_all_files", " --make-bed --out /data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i],"/all_chr",".tag"))




