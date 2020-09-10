#merge the r2 results of AUC
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01)

total <- 4*3*4*3
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
eth.vec <- rep(0,total)
r2.result <- rep(0,total)
l_vec <- rep(0,total)
m_vec <- rep(0,total)
temp = 1
method <- c("eurcoef","tarcoef","eb")
method_vec <- rep("c",total)
for(q in 1:3){
  for(i in 2:5){
    for(l in 1:3){
      for(m in 1:4){
        load(paste0(cur.dir,eth[i],"/r2_eursnp_",method[q],"_rho_",l,"_size_",m))
        eth.vec[temp] = eth[i]
        r2.result[temp] <- mean(r2.vec)
        l_vec[temp] <- l
        m_vec[temp] <- m
        method_vec[temp] <- method[q]
        temp = temp+1
      }
    }
  }
  
}
#best r2 result by varying the p-value threshold
LD.clump.result <- data.frame(eth.vec,r2.result,l_vec,m_vec,method_vec)
colnames(LD.clump.result) <- "r2.vec"
save(LD.clump.result,file = paste0(cur.dir,"best_eur_snp_result.rdata"))
