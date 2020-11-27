#use plink to subset training data for clumping
#load the p-value results
args = commandArgs(trailingOnly = T)
#i represent ethnic group
i = as.numeric(args[[1]])

library(data.table)
library(dplyr)
 sid<-Sys.getenv('SLURM_JOB_ID')
 dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
 temp.dir = paste0('/lscratch/',sid,'/test')
 
 
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
fam <- as.data.frame(fread(paste0(cur.dir,eth[i],"/all_chr.tag.fam")))
ref_fam = fam[1:3000,]
write.table(ref_fam,file = paste0(temp.dir,"/ref_fam.fam"),row.names=F,col.names=F,quote=F)



res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",cur.dir,eth[i],"/all_chr.tag --keep ",temp.dir,"/ref_fam.fam --out ",temp.dir,"/clump_ref_all_chr --make-bed"))
if(res==2){
  stop()
}
system(paste0("mv ",temp.dir,"/clump_ref_all_chr.fam ",cur.dir,eth[i],"/"))
system(paste0("mv ",temp.dir,"/clump_ref_all_chr.bim ",cur.dir,eth[i],"/"))
system(paste0("mv ",temp.dir,"/clump_ref_all_chr.bed ",cur.dir,eth[i],"/"))
#       temp = temp +1 
#     }
#   }
#   
#   
# }
# write.table(code,file = "/data/zhangh24/multi_ethnic/code/LD_simulation_large/LD_clumping.sh",row.names=F,col.names=F,quote=F)
#the p-value takes the min(AFR,EUR)
#if the SNP only exist







#data <- fread2("/gpfs/gsfs11/users/zhangh24/multi_ethnic/result/LD_simulation_new/EUR/temp/chr_13.txt",header=T)








