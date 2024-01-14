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
  for(i in 2:5){
    for(l in 1:3){
      for(m in 1:4){
        
        result.data.temp.list = list()
        for(i_rep in 1:10){
          load( paste0(out.dir,eth[i],
                       "/weighted_prs_five_ci_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".rdata"))
          result.data.temp.list[[i_rep]] = weightedprs.result
        }
        
        result.data.temp = rbindlist(result.data.temp.list)
        
        result.data = as.vector(colMeans(result.data.temp))
        r2_result = data.frame(r2 = result.data[1],
                               r2_low = result.data[2],
                               r2_high = result.data[3])
        eth.vec[temp] = eth[i]
        
        l_vec[temp] <- l
        m_vec[temp] <- m
        ga_vec[temp] <- i1
        result.table.list[[temp]] = data.frame(eth_vec = eth[i],
                                          r2_result,
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
save(result.table,file = paste0(out.dir,"weightedprs.CT.95CI.rdata"))







#merge the r2 results of AUC

eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)

total <- 5*3*4*10
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
eth.vec <- rep(0,total)
r2.vec <- rep(0,total)
l_vec <- rep(0,total)
m_vec = rep(0,total)
ga_vec = rep(0,total)
temp = 1
result.table.list = list()
for(i1 in 1:5){
  for(i in 2:5){
    for(l in 1:3){
      for(m in 1:4){
        
        result.data.temp.list = list()
        for(i_rep in 1:10){
          load( paste0(out.dir,eth[i],
                       "/weighted_prs_five_ci_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".rdata"))
          result.data.temp.list[[i_rep]] = weightedprs.result
        }
        
        result.data.temp = rbindlist(result.data.temp.list)
        
        result.data = result.data.temp[,1,drop=F]
        
        eth.vec = rep(eth[i],n.rep)
        
        l_vec<- rep(l,n.rep)
        m_vec <- rep(m,n.rep)
        ga_vec <- rep(i1,n.rep)
        result.table.list[[temp]] = data.frame(eth_vec = eth.vec,
                                               r2.vec = result.data,
                                               l_vec = l_vec,
                                               m_vec = m_vec,
                                               ga_vec = ga_vec,
                                               rep_vec = c(1:n.rep),
                                               method_vec = rep("Weighted PRS (CT)",n.rep))
        temp = temp+1
      }
    }
    
  }
  
  
}
library(data.table)

result.table = rbindlist(result.table.list)

#LD.clump.result.p.rep)
save(result.table,file = paste0(out.dir,"weightedprs.CT.rep.rdata"))





