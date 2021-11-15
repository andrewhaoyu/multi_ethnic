#merge the by chr prscsx_score
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
eth <- c("EUR","AFR","AMR","EAS","SAS")

library(data.table)




for(i in 3:5){
  setwd(paste0(out.dir.sum,eth[i],"/prscsx/"))
  for(l in 1:3){
    for(m in 1:4){
      for(i_rep in 1:10){
        for(i1 in 1:5){
          result.list = list()
          
          for(j in 1:22){
file = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",eth[i],
              "_pst_eff_a1_b0.5_phi1e-02_chr",j,".txt")
data = fread(file)
result.list[[j]] = data
          }
          result = rbindlist(result.list)
          fwrite(result,file = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",eth[i],
                                           "_pst_eff_a1_b0.5_phi1e-02.txt"),row.names = F,col.names = T,sep = " ")
        }
      }
    }
  }
  system(paste0("rm ",out.dir.sum,eth[i],"/prscsx/*chr*"))
}
