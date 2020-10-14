#use GCTA to generate phenotypes
#generate phenotypes given different genetic architecture
args = commandArgs(trailingOnly = T)
#i represent ethnic group
i = as.numeric(args[[1]])
#l represent causal snps proportion
l = as.numeric(args[[2]])
#i1 represent genetic architecture
i1 = as.numeric(args[[3]])
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"


sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
eth <- c("EUR","AFR","AMR","EAS","SAS")
system(paste0("cp ",cur.dir,eth[i],"/select.cau.snp.bed /lscratch/",sid,"/test/",eth[i],"_select.cau.snp.bed"))
system(paste0("cp ",cur.dir,eth[i],"/select.cau.snp.bim /lscratch/",sid,"/test/",eth[i],"_select.cau.snp.bim"))
system(paste0("cp ",cur.dir,eth[i],"/select.cau.snp.fam /lscratch/",sid,"/test/",eth[i],"_select.cau.snp.fam"))

select.cau <- read.table(paste0(out.dir,eth[i],"/select.cau_rho",l,"_",i1),header=F)

herit <- nrow(select.cau)*var(select.cau$V2)

system(paste0("/data/zhangh24/software/gcta_1.93.2beta/gcta64 --bfile /lscratch/",sid,"/test/",eth[i],"_select.cau.snp --simu-qt --simu-causal-loci ",out.dir,eth[i],"/select.cau_rho",l,"_",i1," --simu-hsq ",herit," --simu-rep 100 --out ",out.dir,eth[i],"/phenotypes_rho",l,"_",i1))
library(data.table)    
all.pheno <- as.data.frame(fread(paste0(out.dir,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
n.train <- c(15000,45000,80000,100000)
fam <- as.data.frame(fread(paste0(cur.dir,eth[i],"/all_chr.tag.fam")))
n <- nrow(fam)

#for(i_rep in cut.point[k]:(cut.point[k+1]-1)){
#take the first ten replicates for plink to run analysis
for(i_rep in 1:10){ 
  print(i_rep)
  pheno <- fam[,1:2]
  for(m in 1:length(n.train)){
    #only use the first set of all pheno as training
    #all the others are used for testing
    
    temp <- all.pheno[,2+i_rep]
    temp[(n.train[m]+1):n] <- NA
    pheno <- cbind(pheno,temp)
  }
  write.table(pheno,file = paste0(out.dir,eth[i],"/pheno_plink_rho_",l,"_rep_",i_rep,"_GA_",i1),row.names = F,col.names = F,quote=F)
  
}









# total <- 5*3*4
# eth_vec <- rep("c",total)
# l_vec <- rep(0,total)
# m_vec <- rep(0,total)
# herit_vec <- rep(0,total)
# eth <- c("EUR","AFR","AMR","EAS","SAS")
# cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
# out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
# eth <- c("EUR","AFR","AMR","EAS","SAS")
# temp <- 1
# 
# for(i in 1:5){
#   for(l in 1:3){
#     #for(m in 1:4){
#       select.cau <- read.table(paste0(out.dir,eth[i],"/select.cau_rho",l,"_",i1),header=F)
#       eth_vec[temp] <- eth[i]
#       l_vec[temp] <- l
#       m_vec[temp] <- m
#       herit_vec[temp] <- nrow(select.cau)*var(select.cau$V2)
#       temp = temp + 1
#     #}
#   }
# }
# herita.table <- data.frame(eth_vec,l_vec,m_vec,herit_vec)
# save(herita.table,file = paste0(out.dir,"herit_table.rdata"))

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
# }s