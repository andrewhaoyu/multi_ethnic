cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
tar.dir <- "/data/DCEG_shared/datashare/simulated_data_meta/summary_data/"


#five ethnics i
eth <- c("EUR","AFR","AMR","EAS","SAS")
for(i in 1:5){
  temp.code <- paste0("cp ",cur.dir,eth[i],"/summary_out_rho_* ",tar.dir,"",eth[i],"/")
  system(temp.code)
}

tar.dir <- "/data/DCEG_shared/datashare/simulated_data_meta/cau_snp_effect_size/"
eth <- c("EUR","AFR","AMR","EAS","SAS")
for(i in 1:5){
  temp.code <- paste0("cp ",cur.dir,eth[i],"/select.cau_rho* ",tar.dir,"",eth[i],"/")
  system(temp.code)
}
