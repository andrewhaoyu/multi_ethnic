#merge the r2 results of AUC
total <- 2*1*3*4
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)


out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
eth.vec <- rep(0,total)
r2.vec <- rep(0,total)
l_vec <- rep(0,total)
m_vec <- rep(0,total)
method_vec <- rep("c",total)
ga_vec <- rep(0,total)
temp = 1
n.rep=n_rep = 3
for(i in 2:5){
  filedir <- paste0(out.dir,eth[i])
  files <- dir(path = filedir,pattern=paste0("r2.list_rho_eb_test_*"),full.names = T)
  for(i1 in 1:2){
    
    #r2.mat <- matrix(0,length(pthres),total)
    
    
    
    for(l in 1:3){
      for(m in 1:1){
        r2.temp <- rep(0,n_rep)
        r2.stack.temp = rep(0,n_rep)
        
        
        
        for(i_rep in 1:n.rep){
          filename = paste0(out.dir,eth[i],"/r2.list_rho_eb_test_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)
          
          load(filename)
          r2.stack.temp[i_rep] = r2.list[[1]]
          
          r2.temp[i_rep] = r2.list[[2]]  
          
          
        }
        eth.vec[temp:(temp+1)] = rep(eth[i],2)
        r2.vec[temp:(temp+1)] <- c(mean(r2.temp),mean(r2.stack.temp))
        l_vec[temp:(temp+1)] <- rep(l,2)
        m_vec[temp:(temp+1)] <- rep(m,2)
        ga_vec[temp:(temp+1)] <- rep(i1,2)
        method_vec[temp:(temp+1)] <- c("TDLD-EB-TEST","TDLD-SLEB-TEST")
        temp = temp+2
      }
    }
  }
  #best r2 result by varying the p-value threshold
  
}  
EB.result.Test  <- data.frame(eth.vec,r2.vec,l_vec,m_vec,method_vec,ga_vec)

save(EB.result.Test,file = paste0(out.dir,"LD.clump.result.EB.test.rdata"))












#merge the r2 results of AUC
total <- 2*1*3*4
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)


out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
eth.vec <- rep(0,total)
r2.vec <- rep(0,total)
l_vec <- rep(0,total)
m_vec <- rep(0,total)
method_vec <- rep("c",total)
ga_vec <- rep(0,total)
temp = 1
n.rep=n_rep = 3
for(i in 2:5){
  filedir <- paste0(out.dir,eth[i])
 
  for(i1 in 1:2){
    
    #r2.mat <- matrix(0,length(pthres),total)
    
    
    
    for(l in 1:3){
      for(m in 1:1){
        r2.temp <- rep(0,n_rep)
        r2.stack.temp = rep(0,n_rep)
        
        
        
        for(i_rep in 1:n.rep){
          filename = paste0(out.dir,eth[i],"/r2.list_rho_eb_test_sd_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)
          
          load(filename)
          r2.stack.temp[i_rep] = r2.list[[1]]
          
          r2.temp[i_rep] = r2.list[[2]]  
          
          
        }
        eth.vec[temp:(temp+1)] = rep(eth[i],2)
        r2.vec[temp:(temp+1)] <- c(mean(r2.temp),mean(r2.stack.temp))
        l_vec[temp:(temp+1)] <- rep(l,2)
        m_vec[temp:(temp+1)] <- rep(m,2)
        ga_vec[temp:(temp+1)] <- rep(i1,2)
        method_vec[temp:(temp+1)] <- c("TDLD-EB-TEST-SD","TDLD-SLEB-TEST-SD")
        temp = temp+2
      }
    }
  }
  #best r2 result by varying the p-value threshold
  
}  
EB.result.Test.sd  <- data.frame(eth.vec,r2.vec,l_vec,m_vec,method_vec,ga_vec)

save(EB.result.Test.sd,file = paste0(out.dir,"LD.clump.result.EB.test.sd.rdata"))

# #r2 result for different p-value threshold
# r2.vec <- rep(0,length(pthres)^2*total)
# l_vec <- rep(0,length(pthres)^2*total)
# m_vec <- rep(0,length(pthres)^2*total)
# eth.vec <- rep(0,length(pthres)^2*total)
# pthres.tar <- rep(0,length(pthres)^2*total)
# pthres.eur <- rep(0,length(pthres)^2*total)
# temp = 1
# num = 0
# for(i in 2:5){
#   for(l in 1:3){
#     for(m in 1:4){
#       r2.temp <- matrix(0,n.rep,length(pthres)^2)
#       for(i_rep in 1:n.rep){
#         load(paste0(out.dir,eth[i],"/r2.list_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
#         r2.temp[i_rep,] = r2.list[[2]]
#       }
#       temp = length(pthres)^2
#       eth.vec[num+(1:temp)] = eth[i]
#       r2.vec[num+(1:temp)] <- colMeans(r2.temp)
#       l_vec[num+(1:temp)] <- l
#       m_vec[num+(1:temp)] <- m
#       pthres.tar[num+(1:temp)] <- expand.grid(pthres,pthres)[,2]
#       pthres.eur[num+(1:temp)] <- expand.grid(pthres,pthres)[,1]
#       num  = num + temp
#     }
#   }
# }
# 
# LD.clump.result.p <- data.frame(eth.vec,r2.vec,l_vec,m_vec,pthres.tar,pthres.eur)
# 
# #cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
# LD.result.list = list(LD.clump.result,LD.clump.result.p)
# 


#}#


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
