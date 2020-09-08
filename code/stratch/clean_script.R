#clean unnecessary files
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
for(i in 2:5){
  for(k1 in 1:12){
    system(paste0("rm ",cur.dir,eth[i],"/prs/prs_two_dim_",k1,"_*.log"))
    system(paste0("rm ",cur.dir,eth[i],"/prs/prs_two_dim_",k1,"_*.nosex"))
    system(paste0("rm ",cur.dir,eth[i],"/prs/prs_two_dim_",k1,"_*.nopred"))
  }
}