#merge the r2 results of AUC
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01)

total <- 5*3*4
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
eth.vec <- rep(0,total)
r2.vec <- rep(0,total)
l_vec <- rep(0,total)
m_vec <- rep(0,total)
r2.mat <- matrix(0,length(pthres),total)
temp = 1
for(i in 1:5){
  for(l in 1:3){
    for(m in 1:4){
      load(paste0(cur.dir,eth[i],"/r2.list_rho_",l,"_size_",m))
      eth.vec[temp] = eth[i]
      r2.vec[temp] <- r2.list[[1]]
      l_vec[temp] <- l
      m_vec[temp] <- m
      
      temp = temp+1
    }
  }
}
#best r2 result by varying the p-value threshold
LD.clump.result <- data.frame(eth.vec,r2.vec,l_vec,m_vec)
#r2 result for different p-value threshold

r2.vec <- rep(0,length(pthres)*total)
l_vec <- rep(0,length(pthres)*total)
m_vec <- rep(0,length(pthres)*total)
eth.vec <- rep(0,length(pthres)*total)
pthres.vec <- rep(0,length(pthres)*total)
temp = 1
num = 0
for(i in 1:5){
  for(l in 1:3){
    for(m in 1:4){
      load(paste0(cur.dir,eth[i],"/r2.list_rho_",l,"_size_",m))
      temp = length(r2.list[[2]])
      eth.vec[num+(1:temp)] = eth[i]
      r2.vec[num+(1:temp)] <- r2.list[[2]]
      l_vec[num+(1:temp)] <- l
      m_vec[num+(1:temp)] <- m
      pthres.vec[num+(1:temp)] <- pthres
      num  = num + temp
    }
  }
}

LD.clump.result.p <- data.frame(eth.vec,r2.vec,l_vec,m_vec,pthres.vec)

cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
LD.result.list = list(LD.clump.result,LD.clump.result.p)
save(LD.result.list,file = paste0(cur.dir,"LD.clump.result.rdata"))