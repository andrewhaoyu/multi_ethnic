#
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/")
#system(paste0("cd /data/zhangh24/multi_ethnic/result/LD_simulation_new/",eth[i1],"/"))
#system(paste0("rm -rf summary_out_rep"))
# system(paste0("rm -rf ",cur.dir,"summary_out_GA"))
# system(paste0("mkdir ",cur.dir,"summary_out_GA"))
# system(paste0("mkdir ",cur.dir,"pheno_GA"))
# system(paste0("cp ",cur.dir,"*.phen ",cur.dir,"pheno_GA/"))
system(paste0("rm -rf ",cur.dir,"pheno_summary_out_GA"))
system(paste0("rm ",cur.dir,"pheno_summary_out_GA.tar.gz"))
system(paste0("mkdir ",cur.dir,"pheno_summary_out_GA"))
for(i1 in 3:5){
  system(paste0("cd ",cur.dir," ;cp phenotypes_rho*_",i1,".phen pheno_summary_out_GA/"))
  system(paste0("cd ",cur.dir," ;cp summary_out_rho*_GA_",i1," pheno_summary_out_GA/"))
}
#system(paste0("cd ",cur.dir," ;mv pheno_GA/* pheno_summary_out_GA/"))
#system(paste0("cd ",cur.dir," ;mv summary_out_GA/* pheno_summary_out_GA/"))
res = system(paste0("cd ",cur.dir," ; tar -zcvf pheno_summary_out_GA.tar.gz pheno_summary_out_GA/"))
if(res==2){
  stop()
}
