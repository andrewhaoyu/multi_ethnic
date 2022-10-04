eth_vec <- c("EUR","AFR","AMR")
trait_vec <-c("height","bmi")
library(data.table)
pthres <- c(Inf,1E-10,5E-08,5E-05,1.0)
result_list = list()
temp = 1
for(l in 1:2){
  for(i in 2:3){
    eth = eth_vec[i]
    trait = trait_vec[l]
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/AOU/clumping_result/PT/",eth,"/",trait)
    file_dir = out.dir
    files = dir(file_dir, pattern  = "CTSLEB_all_test_",full.names = T)
    for(z_ind in 1:length(pthres)){
      eth = eth_vec[i]
      trait = trait_vec[l]
      file_name = paste0(out.dir, "/CTSLEB_all_test_",z_ind,".result")
      if(file_name%in%files){
        load(file_name)
        r2.result = data.frame(eth = eth_vec[i],
                               trait = trait_vec[l],
                               method = paste0("CT-SLEB (pthres =",pthres[z_ind],")" ),
                               r2 = r2_ctsleb
        )
        result_list[[temp]] = r2.result
        temp = temp+ 1  
      }
    }
  }
}

final_result = rbindlist(result_list)
save(final_result, file = "/data/zhangh24/multi_ethnic/result/AOU/CTSLEB/ct_sleb_five_cut_test2.rdata")

#final_result = read.csv("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/PT.csv")
