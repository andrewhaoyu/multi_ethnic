eth = c("EUR", "AFR", "AMR", "EAS", "SAS")
for(i in 1:5){
  system(paste0("cd /data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/ ;",
                "rm -rf prs; mkdir prs"))
}
