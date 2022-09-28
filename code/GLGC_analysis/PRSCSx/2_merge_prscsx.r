
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)

eth <- c("EUR","AFR","AMR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")

for(i in 2:5){
  for(l in 1:4){
    trait = trait_vec[l]
    out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/GLGC/prs/PRSCSX/",eth[i],"/",trait,"/")
    phi = c("1e+00","1e-02","1e-04","1e-06")
    
    for(v in 1:4){
      data.list = list()
      for(j in 1:22){
        data = fread(paste0(out.dir.prs,"sum_",eth[i],"_pst_eff_a1_b0.5_phi",phi[v],"_chr",j,".txt"))
        data.list[[j]] = data
      }  
      data = rbindlist(data.list)
      save(data,file = paste0(out.dir.prs,"sum_",eth[i],"_pst_eff_a1_b0.5_phi",phi[v],".rdata"))
      data.list = list()
      for(j in 1:22){
        data = fread(paste0(out.dir.prs,"sum_EUR_pst_eff_a1_b0.5_phi",phi[v],"_chr",j,".txt"))
        data.list[[j]] = data
      }  
      data = rbindlist(data.list)
      save(data,file = paste0(out.dir.prs,"sum_EUR_pst_eff_a1_b0.5_phi",phi[v],".rdata"))
    }
    
  }
}
