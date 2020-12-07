#merge the r2 results of AUC
for(i1 in 1:2){
  eth <- c("EUR","AFR","AMR","EAS","SAS")
  pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
  n_rep = n.rep = 10
  total <- 5*3*1*3
  out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
  eth.vec <- rep(0,total)
  r2.vec <- rep(0,total)
  #r2.mat <- matrix(0,length(pthres),total)
  #r2.mat.stack <-   matrix(0,length(pthres),total)
  #r2.mat.max <- matrix(0,length(pthres),total)
  method_vec <- rep("c",total)
  temp = 1
  for(i in 1:5){
    for(l in 1:3){
      for(m in 1:1){
        r2.temp <- rep(0,n_rep)
        r2.stack.temp = rep(0,n_rep)
        r2.max.temp = rep(0,n_rep)
        for(i_rep in 1:n_rep){
          load(paste0(out.dir,eth[i],"/r2.list_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
          r2.stack.temp[i_rep] = r2.list[[1]]
          r2.max.temp[i_rep] = r2.list[[2]]
          r2.temp[i_rep] = r2.list[[3]]
        }
        eth.vec[temp:(temp+2)] = rep(eth[i],3)
        r2.vec[temp:(temp+2)] <- c(mean(r2.temp),mean(r2.max.temp),mean(r2.stack.temp))
        l_vec[temp:(temp+2)] <- rep(l,3)
        m_vec[temp:(temp+2)] <- rep(m,3)
        method_vec[temp:(temp+2)] <- c("C+T","C+T max","C+T SL")
        temp = temp+3
      }
    }
    
  }
  #best r2 result by varying the p-value threshold
  LD.clump.result <- data.frame(eth.vec,r2.vec,l_vec,m_vec,method_vec)
  LD.result.list = list(LD.clump.result)
  #LD.clump.result.p.rep)
  save(LD.result.list,file = paste0(out.dir,"LD.clump.result.GA_",i1,".rdata"))
  
}
  
  
  
  
  #r2 result for different p-value threshold
  # r2.vec <- rep(0,length(pthres)*total)
  # l_vec <- rep(0,length(pthres)*total)
  # m_vec <- rep(0,length(pthres)*total)
  # eth.vec <- rep(0,length(pthres)*total)
  # pthres.vec <- rep(0,length(pthres)*total)
  # temp = 1
  # num = 0
  # for(i in 1:5){
  #   for(l in 1:3){
  #     for(m in 1:4){
  #       r2.temp <- matrix(0,n_rep,length(pthres))
  #       for(i_rep in 1:n_rep){
  #         load(paste0(cur.dir,eth[i],"/r2.list_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
  #         r2.temp[i_rep,] = r2.list[[2]]
  #       }
  #       temp = length(pthres)
  #       eth.vec[num+(1:temp)] = eth[i]
  #       r2.vec[num+(1:temp)] <- colMeans(r2.temp)
  #       l_vec[num+(1:temp)] <- l
  #       m_vec[num+(1:temp)] <- m
  #       pthres.vec[num+(1:temp)] <- pthres
  #       num  = num + temp
  #     }
  #   }
  # }
  # 
  # LD.clump.result.p <- data.frame(eth.vec,r2.vec,l_vec,m_vec,pthres.vec)
  # 
  # cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
  # LD.result.list = list(LD.clump.result,LD.clump.result.p)
  # 
 
#}
