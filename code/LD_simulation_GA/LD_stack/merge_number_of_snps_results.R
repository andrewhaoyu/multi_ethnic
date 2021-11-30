out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
library(data.table)
result.list = list()
for(i in 2:5){
  for(l in 1:3){
    for(m in 1:4){
      for(i1 in 1:5){
        load(paste0(out.dir,"n_snp_result_eth_",i,"_rho_",l,"_size_",m,"_ga_",i1,".rdata"))        
      }
    }
  }
}
