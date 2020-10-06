#clean unnecessary files
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"



for(i in 2:5){
  system(paste0("rm ",cur.dir,eth[i],"/*.assoc.linear"))
    system(paste0("rm ",cur.dir,eth[i],"/prs/*.log"))
    system(paste0("rm ",cur.dir,eth[i],"/prs/*.nosex"))
    system(paste0("rm ",cur.dir,eth[i],"/prs/*.nopred"))
 
}


cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation/"
for(i in 1:5){
  system(paste0("rm ",cur.dir,eth[i],"/*.assoc.linear"))
  system(paste0("rm ",cur.dir,eth[i],"/*.log"))
  system(paste0("rm ",cur.dir,eth[i],"/*.nosex"))
  system(paste0("rm ",cur.dir,eth[i],"/prs/*.log"))
  system(paste0("rm ",cur.dir,eth[i],"/prs/*.nosex"))
  system(paste0("rm ",cur.dir,eth[i],"/prs/*.nopred"))
  
}
