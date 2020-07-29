#use GCTA to generate phenotypes
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
system(paste0("/data/zhangh24/software/gcta_1.93.2beta/gcta64 --bfile ",cur.dir,eth[i],"/all_chr.tag --simu-qt --simu-causal-loci ",cur.dir,eth[i],"/select.cau_rho",l, " --simu-hsq 0.4 --out ",cur.dir,eth[i],"/phenotypes_rho",l))











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