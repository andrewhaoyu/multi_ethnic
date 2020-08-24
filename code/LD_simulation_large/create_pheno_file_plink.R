#create eht phenotypes file for plink to run
#i for ethnic group
#l for causal proportion
#k for replicates
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
k = as.numeric(args[[3]])
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"

cut.point <- seq(1,101,20)

library(data.table)    
all.pheno <- as.data.frame(fread(paste0(cur.dir,eth[i],"/phenotypes_rho",l,".phen")))
n.train <- c(15000,45000,80000,100000)
fam <- as.data.frame(fread(paste0(cur.dir,eth[i],"/all_chr.tag.fam")))
n <- nrow(fam)
pheno <- fam[,1:2]
for(i_rep in cut.point[k]:(cut.point[k+1]-1)){
  print(i_rep)
  for(m in 1:length(n.train)){
    #only use the first set of all pheno as training
    #all the others are used for testing
    
    temp <- all.pheno[,2+i_rep]
    temp[(n.train[m]+1):n] <- NA
    pheno <- cbind(pheno,temp)
  }
  write.table(pheno,file = paste0(cur.dir,eth[i],"/pheno_plink_rho_",i_rep,"_",l),row.names = F,col.names = F,quote=F)
  
}

