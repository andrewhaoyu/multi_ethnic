out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
eth <- c("EUR","AFR","AMR","EAS","SAS")


total = 5*3*4*5
eth.vec = rep("c",total)
r2.vec = rep(0,total)
l_vec = rep(0,total)
m_vec = rep(0,total)
ga_vec = rep(0,total)
n.rep = 10
method_vec = rep("SCT",total)
temp = 1
for(i in 1:5){
  for(l in 1:3){
    for(m in 1:4){
      for(i1 in 1:5){
        r2.vec.temp = rep(0,n.rep)
        for(i_rep in 1:10){
          load(paste0(out.dir,eth[i],"/r2.SCT_rho_",l,"_size_",m,"_GA_",i1,"_rep_",i_rep))
          r2.vec.temp[i_rep] = r2
        }
        eth.vec[temp] = eth[i]
        l_vec[temp] = l
        m_vec[temp] = m
        ga_vec[temp] = i1
        r2.vec[temp] = mean(r2.vec.temp)
        temp = temp+1
      }
    }
  }
}
SCT.clump.result <- data.frame(eth.vec,r2.vec,l_vec,m_vec,
                              ga_vec,method_vec)
save(SCT.clump.result,file = "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/LD.clump.result.SCT.rdata")
