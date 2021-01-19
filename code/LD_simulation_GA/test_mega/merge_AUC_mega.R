#merge the r2 results of AUC
for(q in 1:2){
  snplist_vec = c("mega","hm")
  snplist <- snplist_vec[q]
  
  for(i1 in 1:2){
    eth <- c("EUR","AFR","AMR","EAS","SAS")
    pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
    n_rep = 3
    total <- 5*3*4 
    cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
    eth.vec <- rep(0,total)
    r2.vec <- rep(0,total)
    l_vec <- rep(0,total)
    m_vec <- rep(0,total)
    r2.mat <- matrix(0,length(pthres),total)
    
    
    
    temp = 1
    for(i in 1:5){
      allfiles <- dir(path = paste0(cur.dir,eth[i]),pattern = "r2.list_rho_",full.names = T)
      for(l in 1:3){
        for(m in 1:4){
          r2.temp <- rep(0,n_rep)
          for(i_rep in 1:n_rep){
            file <- paste0(cur.dir,eth[i],"/r2.list_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",snplist)
            if(file %in% allfiles==T){
              load(file)
              r2.temp[i_rep] = r2.list[[1]]  
            }else{
              r2.temp[i_rep] = NA
            }
            
          }
          eth.vec[temp] = eth[i]
          r2.vec[temp] <- mean(r2.temp,na.rm = T)
          l_vec[temp] <- l
          m_vec[temp] <- m
          
          temp = temp+1  
        }
      }
    }
    #best r2 result by varying the p-value threshold
    LD.clump.result <- data.frame(eth.vec,r2.vec,l_vec,m_vec)
    
    
    
    
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
    #         load(paste0(cur.dir,eth[i],"/r2.list_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",snplist))
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
    cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
    LD.result.list = list(LD.clump.result)
    #,
                          #LD.clump.result.p)
    
    #LD.result.list = list(LD.clump.result,LD.clump.result.p)
    #LD.clump.result.p.rep)
    save(LD.result.list,file = paste0(cur.dir,"LD.clump.result.GA_",i1,"_",snplist,".rdata"))
  }
  
}











#merge the r2 results of AUC
for(q in 1:1){
  snplist_vec = c("mega","hm")
  snplist <- snplist_vec[q]
  
  for(i1 in 1:2){
    eth <- c("EUR","AFR","AMR","EAS","SAS")
    n_rep = 10
    total <- 5*3*4*n_rep
    cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
    eth.vec <- rep(0,total)
    r2.vec <- rep(0,total)
    l_vec <- rep(0,total)
    m_vec <- rep(0,total)
    rep_vec <- rep(0,total)
    n_rep = 3
    total <- 5*3*4 
    cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
    eth.vec <- rep(0,total)
    r2.vec <- rep(0,total)
    l_vec <- rep(0,total)
    m_vec <- rep(0,total)
    r2.mat <- matrix(0,length(pthres),total)
    
    
    
    temp = 1
    for(i in 1:5){
      allfiles <- dir(path = paste0(cur.dir,eth[i]),pattern = "r2.list_rho_",full.names = T)
      for(l in 1:3){
        for(m in 1:4){
          r2.temp <- rep(0,n_rep)
          for(i_rep in 1:n_rep){
            file <- paste0(cur.dir,eth[i],"/r2.list_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",snplist)
            if(file %in% allfiles==T){
              load(file)
              r2.temp[i_rep] = r2.list[[1]]  
            }else{
              r2.temp[i_rep] = NA
            }
            
          }
          eth.vec[temp] = eth[i]
          r2.vec[temp] <- mean(r2.temp,na.rm = T)
          l_vec[temp] <- l
          m_vec[temp] <- m
          
          temp = temp+1  
        }
      }
    }
    #best r2 result by varying the p-value threshold
    LD.clump.result <- data.frame(eth.vec,r2.vec,l_vec,m_vec)
    
    
    
    
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
    #         load(paste0(cur.dir,eth[i],"/r2.list_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_",snplist))
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
    cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
    LD.result.list = list(LD.clump.result)
    #,
    #LD.clump.result.p)
    
    #LD.result.list = list(LD.clump.result,LD.clump.result.p)
    #LD.clump.result.p.rep)
    save(LD.result.list,file = paste0(cur.dir,"LD.clump.result.GA_",i1,"_",snplist,".rdata"))
  }
  
}
