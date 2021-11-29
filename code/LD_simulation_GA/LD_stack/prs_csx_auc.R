#calculate the AUC for EUR SNPs with EUR coefficients
#with target ethnic coefficients
#with EB coefficients
#calculate AUC for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
#q for three different methods

#i_rep = as.numeric(args[[5]])


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
eth.vec <- rep("c",total)
r2.result <- rep(0,total)
l_vec <- rep(0,total)
m_vec <- rep(0,total)
ga_vc <- rep(0,total)
method_vec <- rep("PRS-CSx",total)
temp= 1
phi = c("1e-02","1e-04","1e-06")

for(i in 2:5){
  for(i1 in 1:5){
    for(l in 1:3){
      #load the phenotype file
      y <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
      y <- y[,2+(1:n.rep)]
      n <- nrow(y)
      y_test_mat <- y[(100000+1):nrow(y),]
      for(m in 1:4){
        #for(q in 1:3){
          
          print(m)
          
        r2.vad.rep <- rep(0,n.rep)
          
           for(i_rep in 1:n.rep){
            r2.test.rep <- rep(0,3)
            coef.mat = matrix(NA,3,2)
             #find the best phi
            for(v in 1:3){
       
              #load target prs
              filename <- paste0(out.dir.sum,eth[i],"/prscsx/prs_csx_",eth[i],"_rho_",l,"_size_",m,"_GA_",i1,"_phi",phi[v],".sscore")
              prs.temp.tar <- as.data.frame(fread(filename))
              #load eur prs
              filename <- paste0(out.dir.sum,eth[i],"/prscsx/prs_csx_EUR_rho_",l,"_size_",m,"_GA_",i1,"_phi",phi[v],".sscore")
              prs.temp.eur <- as.data.frame(fread(filename))
              prs.score.tar.test <- prs.temp.tar[1:10000,4+i_rep]
              prs.score.eur.test <- prs.temp.eur[1:10000,4+i_rep]
              y.test = y_test_mat[1:10000,i_rep]
            
              model1 <- lm(y.test~prs.score.tar.test+prs.score.eur.test)
              coef.mat[v,] = coefficients(model1)[2:3]
              r2.test.rep[v] <- summary(model1)$r.square
              
            }
            idx <- which.max(r2.test.rep)
            filename <- paste0(out.dir.sum,eth[i],"/prscsx/prs_csx_",eth[i],"_rho_",l,"_size_",m,"_GA_",i1,"_phi",phi[idx],".sscore")
            prs.temp.tar <- as.data.frame(fread(filename))
            #load eur prs
            filename <- paste0(out.dir.sum,eth[i],"/prscsx/prs_csx_EUR_rho_",l,"_size_",m,"_GA_",i1,"_phi",phi[idx],".sscore")
            prs.temp.eur <- as.data.frame(fread(filename))
            prs.score.tar.vad <- prs.temp.tar[10001:20000,4+i_rep]
            prs.score.eur.vad <- prs.temp.eur[10001:20000,4+i_rep]
            prs = cbind(prs.score.tar.vad,prs.score.eur.vad)%*%coef.mat[idx,]
            y.vad = y_test_mat[10001:20000,i_rep]
            model = lm(y.vad~prs)
            r2.vad.rep[i_rep] <- summary(model)$r.square
            #for(k in 1:length(pthres)){
            
            #get the number of
            
           
           
            
            
          }
          
          
          eth.vec[temp] <- eth[i]
          r2.result[temp] <- mean(r2.vad.rep)
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

prscsx.result <- data.frame(eth.vec,
                            r2.vec = r2.result,
                            l_vec,
                            m_vec,
                            ga_vec=ga_vc,
                            method_vec = method_vec)

save(prscsx.result,file = paste0(out.dir,
                                 "prscsx.result.rdata"))




