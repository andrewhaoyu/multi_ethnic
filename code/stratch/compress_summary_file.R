#
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i1],"/")
#system(paste0("cd /data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i1],"/"))
#system(paste0("rm -rf summary_out_rep"))
# system(paste0("rm -rf ",cur.dir,"summary_out_GA"))
# system(paste0("mkdir ",cur.dir,"summary_out_GA"))
# system(paste0("mkdir ",cur.dir,"pheno_GA"))
# system(paste0("cp ",cur.dir,"*.phen ",cur.dir,"pheno_GA/"))
res = system(paste0("cd ",cur.dir," ; tar -zcvf pheno_GA.tar.gz pheno_GA/"))
if(res==2){
  stop()
}
res = system(paste0("cd ",cur.dir," ;tar -zcvf summary_out_GA.tar.gz summary_out_GA/"))
if(res==2){
  stop()
}