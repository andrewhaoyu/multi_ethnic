out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])
i1 = as.numeric(args[[5]])
j = as.numeric(args[[6]])

eth = c("EUR","AFR","AMR","EAS","SAS")
result.list = list()
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
            if(file%in%filename==F){
              
            }
          }
        }
      }
    }
  }
}