total <- 3*4*2
eth <- c("EUR","AFR","AMR","EAS","SAS")



out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
eth.vec <- rep(0,total)
r2.vec <- rep(0,total)
l_vec <- rep(0,total)
m_vec <- rep(0,total)
k_vec = rep(0,total)
method_vec <- rep("c",total)
ga_vec <- rep(0,total)
temp = 1
n.rep=n_rep = 10
for(i in 2:2){
  filedir <- paste0(out.dir,eth[i])
  #files <- dir(path = filedir,pattern=paste0("r2.list_rho_eb_*"),full.names = T)
  for(i1 in 1:1){
    
    #r2.mat <- matrix(0,length(pthres),total)
    
    
    
    for(l in 1:3){
      for(m in 1:4){
        
          r2.temp <- rep(0,n_rep)
          r2.test = rep(0,n_rep)
          r2.vad = rep(0,n_rep)
          
          
          for(i_rep in 1:n.rep){
            filename = paste0(out.dir,eth[i],"/r2_overfit_rho_eb_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)
            
            load(filename)
            r2.test[i_rep] = r2.list[1]
            r2.vad[i_rep] = r2.list[2]
            
            
          }
          eth.vec[temp:(temp+1)] = rep(eth[i],2)
          r2.vec[temp:(temp+1)] <- c(mean(r2.test), mean(r2.vad))
          l_vec[temp:(temp+1)] <- rep(l,2)
          m_vec[temp:(temp+1)] <- rep(m,2)
          
          ga_vec[temp:(temp+1)] <- rep(i1,2)
          method_vec[temp:(temp+1)] <- c("Tuning", "Validation")
          temp = temp+2
        
      }
    }
  }
  
  #best r2 result by varying the p-value threshold
  
}  
EB.result <- data.frame(eth.vec,r2.vec,l_vec,m_vec,method_vec,ga_vec)
save(EB.result,file = paste0(out.dir,"CTSLEB_overfitting.rdata"))
