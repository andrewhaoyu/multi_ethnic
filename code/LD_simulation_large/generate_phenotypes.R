#use GCTA to generate phenotypes
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
system(paste0("/data/zhangh24/software/gcta_1.93.2beta/gcta64 --bfile ",cur.dir,eth[i],"/select.cau.snp --simu-qt --simu-causal-loci ",cur.dir,eth[i],"/select.cau_rho",l, " --simu-hsq 0.4 --simu-rep 100 --out ",cur.dir,eth[i],"/phenotypes_rho",l))



# for(i in 1:5){
#   for(l in 1:3){
#   }
# }


#data <- as.data.frame(fread("/data/zhangh24/multi_ethnic/result/LD_simulation_new/EUR/phenotypes_rho3.phen"))
# var(data[,3])
# var(data[,4])
# var(data[,5])
# var(data[,100])
# beta <- read.table(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i],"/select.cau_rho",l))
# var(beta[,2])*nrow(beta)

# for(i in 2:5){
#   fam <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i],"/all_chr.tag.fam")))
#   fam[,1] <- c(1:nrow(fam))
#   fam[,2] <- c(1:nrow(fam))
#   write.table(fam, file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i],"/all_chr.tag.fam"),row.names = F,quote = F,col.names = F)
# }
# 
# i = 1
# for(j in 1:22){
#   fam <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i],"/chr",j,".tag.fam")))
#   fam[,1] <- c(1:nrow(fam))
#   fam[,2] <- c(1:nrow(fam))
#   fam[,6] <- -9
#   write.table(fam, file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i],"/chr",j,".tag.fam"),row.names = F,quote = F,col.names = F)
# }