args = commandArgs(trailingOnly = T)
l = as.numeric(args[[1]])
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)

eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")

#connect chr data into one
phi = c("1e+00","1e-02","1e-04","1e-06")
out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/PRSCSx/EUR/",trait[l],"/")
for(i in 1:5){
  for(v in 1:4){
    data.list = list()
    for(j in 1:22){
      data = fread(paste0(out.dir.prs,"sum_five_",eth[i],"_pst_eff_a1_b0.5_phi",phi[v],"_chr",j,".txt"))
      data.list[[j]] = data
    }  
    data = rbindlist(data.list)
    save(data,file = paste0(out.dir.prs,"sum_five_",eth[i],"_pst_eff_a1_b0.5_phi",phi[v],".rdata"))
    
  }
}
#   data.list = list()
#   for(j in 1:22){
#     data = fread(paste0(out.dir.prs,"sum_EUR_pst_eff_a1_b0.5_phi",phi[v],"_chr",j,".txt"))
#     data.list[[j]] = data
#   }  
#   data = rbindlist(data.list)
#   save(data,file = paste0(out.dir.prs,"sum_EUR_pst_eff_a1_b0.5_phi",phi[v],".rdata"))
# }
