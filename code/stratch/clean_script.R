#clean unnecessary files
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])

eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
system(paste0("rm -rf ",cur.dir,eth[i],"/*mega*"))
# cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
# system(paste0("rm -rf ",cur.dir,eth[i],"/prs"))
#system(paste0("rm -rf /data/zhangh24/multi_ethnic/result/LD_simulation"))

#for(i in 1:5){
  #system(paste0("rm -r ",cur.dir,eth[i],"/prs"))
  #system(paste0("mkdir ",cur.dir,eth[i],"/prs"))
#   system(paste0("rm -rf ",cur.dir,eth[i],"/prs"))
#   system(paste0("rm ",cur.dir,eth[i],"/*.clumped"))
#     system(paste0("rm ",cur.dir,eth[i],"/r2.list*"))
#    # system(paste0("rm ",cur.dir,eth[i],"/*.clumped"))
#     system(paste0("rm ",cur.dir,eth[i],"/summary_out_*"))
#     
#     system(paste0("rm ",cur.dir,eth[i],"/pheno_summary_stat.tar.gz*"))
#     system(paste0("rm -rf ",cur.dir,eth[i],"/pheno_summary_stat"))
#     system(paste0("rm -rf ",cur.dir,eth[i],"/pheno_plink_*"))
#     system(paste0("rm -rf ",cur.dir,eth[i],"/*.phen"))
#     #system(paste0("rm ",cur.dir,eth[i],"/prs/*.log"))
#     #system(paste0("rm ",cur.dir,eth[i],"/prs/*.nosex"))
#     #system(paste0("rm ",cur.dir,eth[i],"/prs/*.nopred"))
#     #system(paste0("rm ",cur.dir,eth[i],"/summary_out_rho_*"))
#     #system(paste0("rm ",cur.dir,eth[i],"/prs/*.nosex"))
#     #system(paste0("rm ",cur.dir,eth[i],"/prs/*.nopred"))
# 
# }
# 

# cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
# 
# 
# 
# for(i in 1:5){
#   system(paste0("rm ",cur.dir,eth[i],"/summary_out_MAF_*"))
#   
#   #system(paste0("rm ",cur.dir,eth[i],"/prs/*.nosex"))
#   #system(paste0("rm ",cur.dir,eth[i],"/prs/*.nopred"))
#   
# }


# cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation/"
# for(i in 1:5){
#   system(paste0("rm ",cur.dir,eth[i],"/*.assoc.linear"))
#   system(paste0("rm ",cur.dir,eth[i],"/*.log"))
#   system(paste0("rm ",cur.dir,eth[i],"/*.nosex"))
#   system(paste0("rm ",cur.dir,eth[i],"/prs/*.log"))
#   system(paste0("rm ",cur.dir,eth[i],"/prs/*.nosex"))
#   system(paste0("rm ",cur.dir,eth[i],"/prs/*.nopred"))
#   
# }
