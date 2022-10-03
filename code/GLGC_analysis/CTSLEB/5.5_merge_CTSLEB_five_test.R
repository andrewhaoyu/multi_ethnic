eth_vec <- c("EUR","AFR","AMR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
pthres <- c(Inf,1E-10,5E-08,5E-05,1.0)
library(data.table)
result_list = list()
temp = 1
for(l in 1:4){
  for(i in 2:5){
    eth = eth_vec[i]
    trait = trait_vec[l]
    out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/",eth,"/",trait)
    file_dir = out.dir
    files = dir(file_dir, pattern  = "CTSLEB_all_test_",full.names = T)
    for(z_ind in 1:length(pthres)){
      eth = eth_vec[i]
      trait = trait_vec[l]
      file_name = paste0(out.dir, "/CTSLEB_all_test_",z_ind,".result")
      if(file_name%in%files){o
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
save(final_result, file = "/data/zhangh24/multi_ethnic/result/GLGC/CTSLEB/ct_sleb_five_cut_test.rdata")

#final_result = read.csv("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/PT.csv")
