#
new.job <- NULL
error.file <- read.table("/data/zhangh24/multi_ethnic/code/LD_simulation_large/eur_snp_error",header=F)

error.file.temp <- error.file[1:880,]
r.job <- read.table("/data/zhangh24/multi_ethnic/code/LD_simulation_large/calculate_prs_eursnp_eurcoef.sh",header=F)
idx <- which(error.file.temp[,3]=="TIMEOUT"|error.file.temp[,3]=="FALIED")
new.job <- rbind(new.job,r.job[idx,])

error.file.temp <- error.file[881:1760,]
r.job <- read.table("/data/zhangh24/multi_ethnic/code/LD_simulation_large/calculate_prs_eursnp_tarcoef.sh",header=F)
idx <- which(error.file.temp[,3]=="TIMEOUT"|error.file.temp[,3]=="FAILED")
new.job <- rbind(new.job,r.job[idx,])


error.file.temp <- error.file[1761:2640,]
r.job <- read.table("/data/zhangh24/multi_ethnic/code/LD_simulation_large/calculate_prs_eursnp_EBcoef.sh",header=F)
idx <- which(error.file.temp[,3]=="TIMEOUT"|error.file.temp[,3]=="FAILED")
new.job <- rbind(new.job,r.job[idx,])

write.table(new.job,file  = "/data/zhangh24/multi_ethnic/code/stratch/EUR_prs_resubmit.sh",
            row.names = F,col.names = F,quote=F)
