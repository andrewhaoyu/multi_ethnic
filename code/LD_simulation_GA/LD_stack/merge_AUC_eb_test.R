#merge the r2 results of AUC
#total <- 5*4*3*4
total <- 5*4*3*4
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01)


out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
eth.vec <- rep(0,total) 
r2.vec <- rep(0,total)
l_vec <- rep(0,total)
m_vec <- rep(0,total)
method_vec <- rep("c",total)
ga_vec <- rep(0,total)
temp = 1
n.rep=n_rep = 10
for(i in 2:5){
  filedir <- paste0(out.dir,eth[i])
  files <- dir(path = filedir,pattern=paste0("r2.list_rho_ebtest_*"),full.names = T)
  for(i1 in 1:5){
    
    #r2.mat <- matrix(0,length(pthres),total)
    
    
    
    for(l in 1:3){
      for(m in 1:4){
        r2.temp <- rep(0,n_rep)
        r2.stack.temp = rep(0,n_rep)
        
        
        
        for(i_rep in 1:n.rep){
          filename = paste0(out.dir,eth[i],"/r2.list_rho_ebtest_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)
          
          load(filename)
          r2.stack.temp[i_rep] = r2.list[[1]]
          
          #r2.temp[i_rep] = r2.list[[2]]  
          
          
        }
        eth.vec[temp:(temp+1)] = rep(eth[i],1)
        r2.vec[temp:(temp+1)] <- c(mean(r2.stack.temp))
        l_vec[temp:(temp+1)] <- rep(l,1)
        m_vec[temp:(temp+1)] <- rep(m,1)
        ga_vec[temp:(temp+1)] <- rep(i1,1)
        method_vec[temp:(temp+1)] <- c("CT-SLEB (two ancestries)")
        temp = temp+1
      }
    }
  }
  #best r2 result by varying the p-value threshold
  
}  
EB.result <- data.frame(eth.vec,r2.vec,l_vec,m_vec,method_vec,ga_vec)
save(EB.result,file = paste0(out.dir,"LD.clump.result.EBtest.rdata"))

total <- 5*4*3*4
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01)






out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
eth.vec <- rep(0,total) 
r2.vec <- rep(0,total)
l_vec <- rep(0,total)
m_vec <- rep(0,total)
method_vec <- rep("c",total)
ga_vec <- rep(0,total)
temp = 1
filescript = rep("c",total)
for(i in 2:5){
  filedir <- paste0(out.dir,eth[i])
  files <- dir(path = filedir,pattern=paste0("r2.list_rho_ebtest_*"),full.names = T)
  for(i1 in 1:5){
    
    #r2.mat <- matrix(0,length(pthres),total)
    
    
    
    for(l in 1:3){
      for(m in 1:4){
        r2.temp <- rep(0,n_rep)
        r2.stack.temp = rep(0,n_rep)
        
        
        
        for(i_rep in 1:n.rep){
          filename = paste0(out.dir,eth[i],"/r2.list_rho_ebtest_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)
          
         if(filename%in%files==F){
           filelist[temp] = paste0("Rscript /data/zhangh24/multi_ethnic/code/LD_simulation_GA/LD_stack/eb_AUC_test.R ",
                                     i," ",l," ",m," ",i_rep," ",i1)
           temp = temp+1
         }
          
          #r2.temp[i_rep] = r2.list[[2]]  
          
          
        }
        
     
      }
    }
  }
  #best r2 result by varying the p-value threshold
  
}  

filescript = filescript[1:(temp-1)]