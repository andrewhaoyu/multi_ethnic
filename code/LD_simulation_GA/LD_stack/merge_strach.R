#merge the r2 results of AUC
library(dplyr)
for(i1 in 1:2){
  eth <- c("EUR","AFR","AMR","EAS","SAS")
  pthres <- c(5E-08,5E-07,5E-06,5E-05,1E-04,1E-03,1E-02,0.5)
  
  total <- 4*3*3
  out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
  eth.vec <- rep(0,total)
  r2.vec <- rep(0,total)
  l_vec <- rep(0,total)
  m_vec <- rep(0,total)
  method_vec <- rep("c",total)
  n.rep=n_rep = 3
  #r2.mat <- matrix(0,length(pthres),total)
  temp = 1
  for(i in 2:5){
    for(l in 1:3){
      for(m in 1:1){
        r2.temp <- rep(0,n_rep)
        r2.stack.temp = rep(0,n_rep)
        r2.max.temp = rep(0,n_rep)
        for(i_rep in 1:n.rep){
          load(paste0(out.dir,eth[i],"/r2.list_rho_eb_test_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
          r2.stack.temp[i_rep] = r2.list[[1]]
          r2.max.temp[i_rep] = r2.list[[2]]
          r2.temp[i_rep] = 0
        }
        eth.vec[temp:(temp+2)] = rep(eth[i],3)
        r2.vec[temp:(temp+2)] <- c(mean(r2.temp),mean(r2.max.temp),mean(r2.stack.temp))
        l_vec[temp:(temp+2)] <- rep(l,3)
        m_vec[temp:(temp+2)] <- rep(m,3)
        method_vec[temp:(temp+2)] <- c("2DLD-eb","2DLD-max-eb","2DLD-SL-eb-test")
        temp = temp+3
      }
    }
  }
  #best r2 result by varying the p-value threshold
  LD.clump.result <- data.frame(eth.vec,r2.vec,l_vec,m_vec,method_vec)
  LD.clump.result = LD.clump.result %>% 
    filter(method_vec=="2DLD-SL-eb-test")
  LD.result.list = list(LD.clump.result)
  save(LD.result.list,file = paste0(out.dir,"LD.clump.result.eb_test_GA_",i1,".rdata"))
}  

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
