#calculate the AUC for EUR SNPs with EUR coefficients
#with target ethnic coefficients
#with EB coefficients
#calculate AUC for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
#q for three different methods

#i_rep = as.numeric(args[[5]])

method <- c("Weighted-PRS")

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
method_vec <- rep("c",total)
temp= 1
for(i in 2:5){
  for(i1 in 1:5){
    for(l in 1:3){
      #load the phenotype file
      y <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
      y <- y[,2+(1:n.rep)]
      n <- nrow(y)
      y_test_mat <- y[(100000+1):110000,]
      y_vad_mat <- y[(110001):120000,]
      for(m in 1:4){
        print(m)
          
          r2.test.rep <- rep(0,n.rep)
          
          for(i_rep in 1:n.rep){
            
            
            #for(k in 1:length(pthres)){
            
            #get the number of
            filename <- paste0(out.dir,eth[i],"/prs/prs_besteur_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".sscore")
            prs.temp <- as.data.frame(fread(filename))
            #load best EUR PRS
            q = 1
            prs.score.eur<- prs.temp[,4+q]
            #load AUC performance for target
            load(paste0(out.dir,eth[i],"/r2.list_rho_",l,"_size_",m,"_GA_",i1))
            best.auc.result <- as.data.frame(r2.list[[2]]) %>% 
              filter(r2.vec.test==max(r2.vec.test))
            r_ind  = best.auc.result$r2_ind_vec
            w_ind = best.auc.result$wc_ind_vec
            k =which(pthres==best.auc.result$pthres_vec)
            filename <- paste0(out.dir,eth[i],"/prs/prs_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".p_value_",k,".profile")
            prs.temp <- fread(filename)
            prs.score.tar <- prs.temp$SCORE
            
            y.test = y_test_mat[,i_rep]
            y.vad = y_vad_mat[,i_rep]
            prs.score = cbind(prs.score.eur,prs.score.tar)
            prs.score.test = prs.score[1:10000,]
            prs.score.vad = prs.score[10001:20000,]
            #y.vad = y_test_mat[(n.test+1):(nrow(y_test_mat)),i_rep]
            model1 <- lm(y.test~prs.score.test)
            beta = coefficients(model1)[2:3]
            model2 = lm(y.vad~prs.score.vad%*%beta)
            
            
            r2.test.rep[i_rep] <- summary(model2)$r.square
            
            
            
          }
          
          
          eth.vec[temp] <- eth[i]
          r2.result[temp] <- mean(r2.test.rep)
          l_vec[temp] <- l
          m_vec[temp] <- m
          ga_vc[temp] <- i1
          method_vec[temp] <- method
          temp = temp+1
          
        }
      }
    }
    
    
  }


weightedprs.result <- data.frame(eth.vec,
                            r2.vec = r2.result,
                            l_vec,
                            m_vec,
                            ga_vec=ga_vc,
                            method_vec = method_vec)

save(weightedprs.result,file = paste0(out.dir,
                                 "weightedprs.result.rdata"))







#load(paste0(out.dir,"LD.clump.result.GA_",i1,".rdata"))


save(r2.vec,file = paste0(out.dir,eth[i],"/r2_eursnp_",method[q],"_rho_",l,"_size_",m,"_GA_",i1))
#write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump.auc.csv")
