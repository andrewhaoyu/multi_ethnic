#merge the r2 results of AUC
for(i1 in 1:2){
  eth <- c("EUR","AFR","AMR","EAS","SAS")
  pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
  
  total <- 4*3*4
  out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
  eth.vec <- rep(0,total)
  r2.vec <- rep(0,total)
  l_vec <- rep(0,total)
  m_vec <- rep(0,total)
  n.rep = 3
  #r2.mat <- matrix(0,length(pthres),total)
  temp = 1
  for(i in 2:5){
    for(l in 1:3){
      for(m in 1:4){
        r2.temp <- rep(0,n.rep)
        for(i_rep in 1:n.rep){
          load(paste0(out.dir,eth[i],"/r2.list_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
          r2.temp[i_rep] = r2.list[[1]]
        }
        eth.vec[temp] = eth[i]
        r2.vec[temp] <- mean(r2.temp)
        l_vec[temp] <- l
        m_vec[temp] <- m
        
        temp = temp+1  
      }
    }
  }
  #best r2 result by varying the p-value threshold
  LD.clump.result <- data.frame(eth.vec,r2.vec,l_vec,m_vec)
  
  
  #r2 result for different p-value threshold
  r2.vec <- rep(0,length(pthres)^2*total)
  l_vec <- rep(0,length(pthres)^2*total)
  m_vec <- rep(0,length(pthres)^2*total)
  eth.vec <- rep(0,length(pthres)^2*total)
  pthres.tar <- rep(0,length(pthres)^2*total)
  pthres.eur <- rep(0,length(pthres)^2*total)
  temp = 1
  num = 0
  for(i in 2:5){
    for(l in 1:3){
      for(m in 1:4){
        r2.temp <- matrix(0,n.rep,length(pthres)^2)
        for(i_rep in 1:n.rep){
          load(paste0(out.dir,eth[i],"/r2.list_two_dim_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
          r2.temp[i_rep,] = r2.list[[2]]
        }
        temp = length(pthres)^2
        eth.vec[num+(1:temp)] = eth[i]
        r2.vec[num+(1:temp)] <- colMeans(r2.temp)
        l_vec[num+(1:temp)] <- l
        m_vec[num+(1:temp)] <- m
        pthres.tar[num+(1:temp)] <- expand.grid(pthres,pthres)[,2]
        pthres.eur[num+(1:temp)] <- expand.grid(pthres,pthres)[,1]
        num  = num + temp
      }
    }
  }
  
  LD.clump.result.p <- data.frame(eth.vec,r2.vec,l_vec,m_vec,pthres.tar,pthres.eur)
  
  #cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
  LD.result.list = list(LD.clump.result,LD.clump.result.p)
  
  save(LD.result.list,file = paste0(out.dir,"LD.clump.result.two.dim_GA_",i1,".rdata"))
  
}

  
LD.clump.result.p <- LD.result.list[[2]]
temp = LD.clump.result.p %>% filter(eth.vec=="SAS"&
                                      l_vec ==3&
                                      m_vec==1)
which.max(LD.clump.result.p$r2)
  #r2 result for different p-value threshold and different rep
  # n_rep =10
  # r2.vec <- rep(0,length(pthres)^2*total*n_rep)
  # l_vec <- rep(0,length(pthres)^2*total*n_rep)
  # m_vec <- rep(0,length(pthres)^2*total*n_rep)
  # eth.vec <- rep(0,length(pthres)^2*total*n_rep)
  # pthres.tar <- rep(0,length(pthres)^2*total*n_rep)
  # pthres.eur <- rep(0,length(pthres)^2*total*n_rep)
  # rep.vec <- rep(0,length(pthres)^2*total*n_rep)
  # temp = 1
  # num = 0
  # for(i in 2:5){
  #   for(l in 1:3){
  #     for(m in 1:4){
  #       for(i_rep in 1:10){
  #         
  #         load(paste0(cur.dir,eth[i],"/r2.list_rho_",l,"_size_",m,"_rep_",i_rep))
  #         temp = length(pthres)
  #         eth.vec[num+(1:temp)] = eth[i]
  #         r2.vec[num+(1:temp)] <- r2.list[[2]]
  #         l_vec[num+(1:temp)] <- l
  #         m_vec[num+(1:temp)] <- m
  #         pthres.vec[num+(1:temp)] <- pthres
  #         rep.vec[num+(1:temp)] <- i_rep
  #         num  = num + temp
  #         
  #       }
  #     }
  #   }
  # }
  # 
  # LD.clump.result.p.rep <- data.frame(eth.vec,r2.vec,l_vec,m_vec,pthres.vec,rep.vec)
  # 
  # cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
  # LD.result.list = list(LD.clump.result,LD.clump.result.p,
  #                       LD.clump.result.p.rep)
  