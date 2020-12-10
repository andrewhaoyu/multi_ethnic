#run GWAS using plink2
#use plink2 to run
args  = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#i_rep represent simulation replication
#i1 represent the genetic architecture
i1 = as.numeric(args[[1]])
# i_rep = as.numeric(args[[2]])
l = as.numeric(args[[2]])
i = as.numeric(args[[3]])
j = as.numeric(args[[4]])



#load the phenotpypes data and use plink to run
for(i_rep in 1:10){
  eth <- c("EUR","AFR","AMR","EAS","SAS")
  library(data.table)
  cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
  out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
  
  sid<-Sys.getenv('SLURM_JOB_ID')
  dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
  temp.dir = paste0('/lscratch/',sid,'/test/')
  eth <- c("EUR","AFR","AMR","EAS","SAS")
  system(paste0("cp ", cur.dir,eth[i],"/chr",j,"train.mega.bed /lscratch/",sid,"/test/",eth[i],"_chr",j,".train.mega.bed"))
  system(paste0("cp ", cur.dir,eth[i],"/chr",j,"train.mega.bim /lscratch/",sid,"/test/",eth[i],"_chr",j,".train.mega.bim"))
  system(paste0("cp ", cur.dir,eth[i],"/chr",j,"train.mega.fam /lscratch/",sid,"/test/",eth[i],"_chr",j,".train.mega.fam"))
  
  res <- system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile /lscratch/",sid,"/test/",eth[i],"_chr",j,".train.mega --out ",temp.dir,"summary_chr_",j,"_rho_",l,"_rep_",i_rep,"_GA_",i1,".out --linear --all-pheno --allow-no-sex --pheno ",out.dir,eth[i],"/pheno_plink_rho_",l,"_rep_",i_rep,"_GA_",i1))
  if(res==2){
    stop()
  }
  system(paste0("mv ",temp.dir,"/*.linear ",out.dir,eth[i],"/."))
  system(paste0('rm -rf /lscratch/',sid,'/test'))
  
}
# for(i in 1:5){
#   fam <- data.frame(fread(paste0(cur.dir,eth[i],"/all_chr.tag.fam")))
#     for(j in 1:22){
#     write.table(fam, file = paste0(cur.dir,eth[i],"/chr",j,".tag.fam"),row.names = F,col.names = F,quote=F)
#   }
# }
# 
# 
