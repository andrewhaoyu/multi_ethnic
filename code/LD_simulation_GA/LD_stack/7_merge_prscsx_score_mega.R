#Goal: merge different CHRs into one single file
args = commandArgs(trailingOnly = T)
# v = as.numeric(args[[1]])
# i = as.numeric(args[[2]])
# l = as.numeric(args[[3]])
# m = as.numeric(args[[4]])
v = as.numeric(args[[1]])
i = 2
l = 1
m = 1
#i_rep = as.numeric(args[[1]])
i1 = 1

#merge the by chr prscsx_score
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
eth <- c("EUR","AFR","AMR","EAS","SAS")

library(data.table)


phi = c("1e+00","1e-02","1e-04","1e-06")

# for(v in 1:3){
#   for(i in 2:5){
setwd(paste0(out.dir.sum,eth[i],"/prscsx_mega/"))
# for(l in 1:3){
#   for(m in 1:4){
for(v in 1:4){
  for(i_rep in 1:10){
    #for(i1 in 1:5){
    result.list.tar = list()
    result.list.eur = list()
    print(i_rep)
    for(j in 1:22){
      file = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",eth[i],
                    "_pst_eff_a1_b0.5_phi",phi[v],"_chr",j,".txt")
      data = fread(file)
      result.list.tar[[j]] = data
      file = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_EUR_pst_eff_a1_b0.5_phi",phi[v],"_chr",j,".txt")
      data = fread(file)
      result.list.eur[[j]] = data
    }
    result = rbindlist(result.list.tar)
    fwrite(result,file = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",eth[i],
                                "_pst_eff_a1_b0.5_phi",phi[v],".txt"),row.names = F,col.names = T,sep = " ")
    result = rbindlist(result.list.eur)
    fwrite(result,file = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_EUR_pst_eff_a1_b0.5_phi"
                                ,phi[v],".txt"),row.names = F,col.names = T,sep = " ")
    
    #}
  }
  
}
#       }
#     }
#     #system(paste0("rm ",out.dir.sum,eth[i],"/prscsx/*chr*"))
#   }
#   
# }
