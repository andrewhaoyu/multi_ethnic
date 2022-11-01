
#merge the by chr prscsx_score
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
eth <- c("EUR","AFR","AMR","EAS","SAS")

out.dir = paste0(out.dir.sum,eth[1],"/prscsx")
#rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1
library(data.table)


phi = c("1e+00","1e-02","1e-04","1e-06")
# for(v in 1:3){
#   for(i in 2:5){

files = dir(path = out.dir, pattern = paste0("rho"))
resubmit_code = rep("c",264000)
temp = 1

 for(l in 1:3){
   for(m in 1:3){
for(i_rep in 1:10){
  print(i_rep)
  for(i1 in 1:5){
    for(j in 1:22){
      for(k in 1:4){
        
        #for(i in 1:5){
          file = paste0("update_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",eth[1],
                        "_pst_eff_a1_b0.5_phi",phi[k],"_chr",j,".txt")
          if(file%in%files==F){
            resubmit_code[[temp]] = paste0("Rscript /data/zhangh24/multi_ethnic/code/LD_simulation_GA/LD_stack/PRSCSx/2.5_rerun_prscsx_five.R ",
                                            l, " " , m, " ", i_rep, " ", i1, " ", k, " ",j)
            temp = temp + 1
          }
        #}
      }
    }
  }
}
   }
 }
resubmit_code = resubmit_code[1:(temp-1)]
write.table(resubmit_code, file = "/data/zhangh24/multi_ethnic/code/LD_simulation_GA/LD_stack/PRSCSx/2.5_rerun_prscsx_five.sh", 
            row.names = F, col.names = F, quote = F)
#     result.list.tar = list()
#     result.list.eur = list()
#     
#     for(j in 1:22){
#       file = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",eth[i],
#                     "_pst_eff_a1_b0.5_phi",phi[k],"_chr",j,".txt")
#       data = fread(file)
#       result.list.tar[[j]] = data
#       file = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_EUR_pst_eff_a1_b0.5_phi",phi[v],"_chr",j,".txt")
#       data = fread(file)
#       result.list.eur[[j]] = data
#     }
#     result = rbindlist(result.list.tar)
#     fwrite(result,file = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",eth[i],
#                                 "_pst_eff_a1_b0.5_phi",phi[v],".txt"),row.names = F,col.names = T,sep = " ")
#     result = rbindlist(result.list.eur)
#     fwrite(result,file = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_EUR_pst_eff_a1_b0.5_phi"
#                                 ,phi[v],".txt"),row.names = F,col.names = T,sep = " ")
#     
#   }
# }
#       }
#     }
#     #system(paste0("rm ",out.dir.sum,eth[i],"/prscsx/*chr*"))
#   }
#   
# }
