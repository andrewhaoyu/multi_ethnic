#merge the r2 results of AUC
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)

total <- 4*3*4
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
eth.vec <- rep(0,total)
r2.vec <- rep(0,total)
l_vec <- rep(0,total)
m_vec <- rep(0,total)
#r2.mat <- matrix(0,length(pthres),total)
temp = 1
for(i in 2:5){
  for(l in 1:3){
    for(m in 1:4){
      r2.temp <- rep(0,10)
      for(i_rep in 1:10){
        load(paste0(cur.dir,eth[i],"/r2.list_eb_rho_",l,"_size_",m,"_rep_",i_rep))
        if(length(r2.list[[1]])==0){
          r2.temp[i_rep] = NA  
        }else{
          r2.temp[i_rep] = r2.list[[1]]
        }
        
      }
      eth.vec[temp] = eth[i]
      r2.vec[temp] <- mean(r2.temp,na.rm=T)
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
      r2.temp <- matrix(0,10,length(pthres)^2)
      for(i_rep in 1:10){
        load(paste0(cur.dir,eth[i],"/r2.list_two_dim_rho_",l,"_size_",m,"_rep_",i_rep))
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

cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
LD.result.list = list(LD.clump.result,LD.clump.result.p)




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
#LD.result.list = list(LD.clump.result)
save(LD.result.list,file = paste0(cur.dir,"LD.clump.result.eb.rdata"))
