args = commandArgs(trailingOnly = T)
k = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])

#merge the by chr prscsx_score
eth <- c("EUR","AFR","AMR","EAS","SAS")
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir = paste0(out.dir.sum,eth[1],"/prscsx")

library(data.table)


phi = c("1e+00","1e-02","1e-04","1e-06")
files = dir(path = out.dir, pattern = paste0("rho"))
# for(v in 1:3){
#   for(i in 2:5){
setwd(out.dir)
# for(l in 1:3){
#   for(m in 1:4){
#for(i_rep in 1:10){
  for(i1 in 1:5){
    for(i in 1:5){
      result.list.tar = list()
    for(j in 1:22){
        if(m==4){
          file = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",eth[i],
                        "_pst_eff_a1_b0.5_phi",phi[k],"_chr",j,".txt")
          data = fread(file)
          result.list.tar[[j]] = data
          
        }else{
          # file = paste0("update_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",eth[i],
          #               "_pst_eff_a1_b0.5_phi",phi[k],"_chr",j,".txt")
          # if(file%in%files){
            data = fread(file)
            result.list.tar[[j]] = data
          #}
          
          
          
        }
        
        
      }
      result = rbindlist(result.list.tar)
      if(m==4){
        fwrite(result,file = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",eth[i],
                                    "_pst_eff_a1_b0.5_phi",phi[k],".txt"),row.names = F,col.names = T,sep = " ")
        
      }else{
        fwrite(result,file = paste0("update_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",eth[i],
                                    "_pst_eff_a1_b0.5_phi",phi[k],".txt"),row.names = F,col.names = T,sep = " ")
        
      }
      
    
  }
  }
#}
#       }
#     }
#     #system(paste0("rm ",out.dir.sum,eth[i],"/prscsx/*chr*"))
#   }
#   
# }
