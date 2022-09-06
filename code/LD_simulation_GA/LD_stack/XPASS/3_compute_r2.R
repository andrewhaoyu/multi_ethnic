#compute R2 for XPASS


library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
#n <- 120000
n.train.vec <- c(15000,45000,80000,100000)
n.rep = 10

#r2 mat represent the r2 matrix for the testing dataset
#column represent the ethnic groups
#row represent different p-value threshold

cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
setwd("/data/zhangh24/multi_ethnic/")


total <- 4*3*4*5
#total <- 4*3*4*1
eth.vec <- rep("c",total)
r2.result <- rep(0,total)
l_vec <- rep(0,total)
m_vec <- rep(0,total)
ga_vc <- rep(0,total)
method_vec <- rep("XPASS",total)
temp= 1
#phi = c("1e+00","1e-02","1e-04","1e-06")

mis_vec_list = NULL 
  mis_temp = 1
for(i in 2:5){
  file_dir = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/xpass")
  files = dir(file_dir, pattern = paste0(".sscore"),
              full.names = T)
  for(i1 in 1:1){
    for(l in 1:3){
# for(i in 2:2){
#   for(i1 in 1:1){
#     for(l in 1:1){
      
      
      #load the phenotype file
      y <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
      y <- y[,2+(1:n.rep)]
      n <- nrow(y)
      y_test_mat <- y[(100000+1):nrow(y),]
      #for(m in 1:4){
      for(m in 1:4){
        #for(q in 1:3){
        
        print(m)
        
        r2.vad.rep <- rep(0,n.rep)
        
        for(i_rep in 1:n.rep){
          r2.test.rep <- rep(0,3)
          coef.mat = matrix(NA,4,2)
          #find the best phi
          
            
            #load target prs
            file_name = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/xpass/rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_PRS.sscore")
            if(file_name%in%files){
              PRS <- as.data.frame(fread(paste0(file_name)))
              prs.tar = PRS$SCORE3_SUM[10001:20000]
              y.vad = y_test_mat[10001:20000,i_rep]
              model = lm(y.vad~prs.tar)
              r2.vad.rep[i_rep] <- summary(model)$r.square
              
            }else{
              mis_vec_list[[mis_temp]] = data.frame(i, l, m, i_rep, i1)
              mis_temp = mis_temp + 1
              r2.vad.rep[i_rep] = NA
            }
          #for(k in 1:length(pthres)){
          
          #get the number of
          
          
          
          
          
        }
        
        
        eth.vec[temp] <- eth[i]
        r2.result[temp] <- mean(r2.vad.rep, na.rm = T)
        l_vec[temp] <- l
        m_vec[temp] <- m
        ga_vc[temp] <- i1
        #method_vec[temp] <- method[q]
        temp = temp+1
        
      }
      # }
    }
    
    
  }
}
  mis_vec = rbindlist(mis_vec_list)
xpass.result <- data.frame(eth.vec,
                            r2.vec = r2.result,
                            l_vec,
                            m_vec,
                            ga_vec=ga_vc,
                            method_vec = method_vec)

save(xpass.result,file = paste0(out.dir,
                                 "xpass.result.rdata"))















for(i in 2:5){
  file_dir = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/xpass")
  files = dir(file_dir, pattern = paste0(".sscore"),
              full.names = T)
  for(i1 in 1:5){
    for(l in 1:3){
      for(m in 1:4){
        #for(q in 1:3){
    
        
        for(i_rep in 1:n.rep){
          
          
          #load target prs
          file_name = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/xpass/rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_PRS.sscore")
          if(file_name%in%files){
            
          }else{
            mis_vec_list[[mis_temp]] = data.frame(i, l, m, i_rep, i1)
            mis_temp = mis_temp + 1
          }
          #for(k in 1:length(pthres)){
          
          #get the number of
          
          
          
          
          
        }
        
    
      }
      # }
    }
    
    
  }
}
mis_vec = rbindlist(mis_vec_list)

run_code = rep("c",nrow(mis_vec))

for(k in 1:nrow(mis_vec)){
  run_code[k] = paste0("Rscript /gpfs/gsfs11/users/zhangh24/multi_ethnic/code/LD_simulation_GA/LD_stack/XPASS/2.5_XPASS_PRS_rerun.R ",
                       mis_vec[k,1]," ",
                       mis_vec[k,2]," ",
                       mis_vec[k,3]," ",
                       mis_vec[k,4]," ",
                       mis_vec[k,5])
}
write.table(run_code, file = "/data/zhangh24/multi_ethnic/code/LD_simulation_GA/LD_stack/XPASS/2.5_XPASS_PRS_rerun.sh",row.names = F, col.names = F, quote=F)
