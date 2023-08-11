#merge the r2 results of AUC

eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)

total <- 5*3*4
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
eth.vec <- rep(0,total)
r2.vec <- rep(0,total)
l_vec <- rep(0,total)
m_vec = rep(0,total)
ga_vec = rep(0,total)
temp = 1
result.table.list = list()
for(i1 in 1:5){
  for(i in 1:5){
    for(l in 1:3){
      for(m in 1:4){
        
        
        
        load(paste0(out.dir,eth[i],"/ct95_",l,"_size_",m,"_GA_",i1))
        
        
        eth.vec[temp] = eth[i]
       
        l_vec[temp] <- l
        m_vec[temp] <- m
        ga_vec[temp] <- i1
        result.table.list[[temp]] = data.frame(eth_vec = eth[i],
                                              result.data,
                                               l_vec = l,
                                               ga_vec = i1,
                                               m_vec = m)
        temp = temp+1
      }
    }
    
  }
  
  
}
library(data.table)

result.table = rbindlist(result.table.list)

#LD.clump.result.p.rep)
save(LD.result.list,file = paste0(out.dir,"LD.clump.result.CT.95CI.rdata"))






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
