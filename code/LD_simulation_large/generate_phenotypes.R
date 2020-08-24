#use GCTA to generate phenotypes
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
system(paste0("/data/zhangh24/software/gcta_1.93.2beta/gcta64 --bfile ",cur.dir,eth[i],"/select.cau.snp --simu-qt --simu-causal-loci ",cur.dir,eth[i],"/select.cau_rho",l, " --simu-hsq 0.4 --simu-rep 100 --out ",cur.dir,eth[i],"/phenotypes_rho",l))

library(data.table)    
all.pheno <- as.data.frame(fread(paste0(cur.dir,eth[i],"/phenotypes_rho",l,".phen")))
n.train <- c(15000,45000,80000,100000)
fam <- as.data.frame(fread(paste0(cur.dir,eth[i],"/all_chr.tag.fam")))
n <- nrow(fam)

#for(i_rep in cut.point[k]:(cut.point[k+1]-1)){
for(i_rep in 1:100){ 
  print(i_rep)
  pheno <- fam[,1:2]
  for(m in 1:length(n.train)){
    #only use the first set of all pheno as training
    #all the others are used for testing
    
    temp <- all.pheno[,2+i_rep]
    temp[(n.train[m]+1):n] <- NA
    pheno <- cbind(pheno,temp)
  }
  write.table(pheno,file = paste0(cur.dir,eth[i],"/pheno_plink_rho_",l,"_rep_",i_rep),row.names = F,col.names = F,quote=F)
  
}



# for(i in 1:5){
#   for(l in 1:3){
#     for(i_rep in 1:100){
#       
#       system(paste0("mv ",cur.dir,eth[i],"/pheno_plink_rho_",i_rep,"_",l," ",cur.dir,eth[i],"/pheno_plink_rho_",l,"_rep_",i_rep))    }
#   }
# }
# 

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