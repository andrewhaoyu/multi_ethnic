#create eht phenotypes file for plink to run
#i for ethnic group
#l for causal proportion
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])

eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
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
