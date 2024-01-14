out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
result_list = list()
temp = 1
for(i in 2:5){
  for(i1 in 1:5){
    for(l in 1:3){
      for(m in 1:4){
        load(paste0(out.dir,
                    "prscsx_five.result_rep_ga_",i1,"_eth_",i,"_rho_",l,
                    "_size_",m,".rdata"))
        result_list[[temp]] = prscsx.result
        temp = temp + 1
      }
    }
  }
}
prscsx.result = rbindlist(result_list)
save(prscsx.result, file = paste0(out.dir,
            "prscsx_five.result_rep_ga.rdata"))
