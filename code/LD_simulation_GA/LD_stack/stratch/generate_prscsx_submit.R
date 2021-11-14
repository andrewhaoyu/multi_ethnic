eth = c("EUR","AFR","AMR","EAS","SAS")
result.list = list()
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
temp = 1
for(i in 2:5){
  file.path = paste0(out.dir.sum,eth[i],"/prscsx/")
  file = dir(file.path,pattern="rho")
  
    for(l in 1:3){
    for(m in 1:4){
      for(i_rep in 1:10){
        for(i1 in 1:5){
          for(j in 1:22){
            filename = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",eth[i],
                              "_pst_eff_a1_b0.5_phi1e-02_chr",j,".txt")
            if(filename%in%file==F){
              temp.file = paste0("Rscript /data/zhangh24/multi_ethnic/code/LD_simulation_GA/LD_stack/stratch/PRScsx_run_temp.R ",i," ",l," ",m," ",i_rep," ",i1," ",j)
              result.list[[temp]] = temp.file
              temp = temp + 1
            }
          }
        }
      }
    }
  }
}
library(data.table)
result.file = data.frame()
for(k in 1:(temp-1)){
  result.file[k,1] = result.list[[k]]
}
write.table(result.file,file = "/data/zhangh24/multi_ethnic/code/LD_simulation_GA/LD_stack/stratch/generate_prscsx_submit.sh",
            row.names = F,col.names = F,quote=F)
