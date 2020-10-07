#run GWAS using plink2
#use plink2 to run
args  = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#i_rep represent simulation replication
#i1 represent the genetic architecture
i = as.numeric(args[[1]])
j = as.numeric(args[[2]])
l = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])
i1 = as.numeric(args[[5]])
#load the phenotpypes data and use plink to run
eth <- c("EUR","AFR","AMR","EAS","SAS")
library(data.table)
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"

sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
eth <- c("EUR","AFR","AMR","EAS","SAS")
system(paste0("cp ", cur.dir,eth[i],"/chr",j,".tag.bed /lscratch/",sid,"/test/",eth[i],"_chr",j,".tag.bed"))
system(paste0("cp ", cur.dir,eth[i],"/chr",j,".tag.bim /lscratch/",sid,"/test/",eth[i],"_chr",j,".tag.bim"))
system(paste0("cp ", cur.dir,eth[i],"/chr",j,".tag.fam /lscratch/",sid,"/test/",eth[i],"_chr",j,".tag.fam"))

system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile /lscratch/",sid,"/test/",eth[i],"_chr",j,".tag --out ",out.dir,eth[i],"/summary_chr_",j,"_rho_",l,"_rep_",i_rep,"_GA_",i1,".out --linear --all-pheno --allow-no-sex --pheno ",out.dir,eth[i],"/pheno_plink_rho_",l,"_rep_",i_rep,"_GA_",i1))


# for(i in 1:5){
#   fam <- data.frame(fread(paste0(cur.dir,eth[i],"/all_chr.tag.fam")))
#     for(j in 1:22){
#     write.table(fam, file = paste0(cur.dir,eth[i],"/chr",j,".tag.fam"),row.names = F,col.names = F,quote=F)
#   }
# }
# 
# 
