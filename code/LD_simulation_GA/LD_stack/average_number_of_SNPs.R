eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
old.out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
library(dplyr)
library(data.table)
setwd("/data/zhangh24/multi_ethnic/")

pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)


i = as.numeric(args[[2]])
l = as.numeric(args[[3]])
m = as.numeric(args[[4]])
i1 = as.numeric(args[[5]])
r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)

i_rep = 1
temp = 1
N_vec = rep(0,10000)
for(i in 2:5){
  for(l in 1:3){
    for(m in 1:4){
      for(i1 in 1:5){
        for(r_ind in 1:length(r2_vec)){
          wc_vec = wc_base_vec/r2_vec[r_ind]
          for(w_ind in 1:length(wc_vec)){
        LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
        clump.snp <- LD[,1,drop=F] 
        N_vec[temp] = nrow(clump.snp)
        temp = temp+1
          }
        }
      }
    }
  }
}
N_vec = N_vec[1:(temp-1)]
mean(N_vec)
