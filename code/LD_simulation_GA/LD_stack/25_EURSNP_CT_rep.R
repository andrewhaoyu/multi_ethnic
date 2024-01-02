#calculate the AUC for EUR SNPs with EUR coefficients
#with target ethnic coefficients
#with EB coefficients
#calculate AUC for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
#q for three different methods

#i_rep = as.numeric(args[[5]])

method <- c("Best EUR SNP (CT)")

library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
#n <- 120000
n.train.vec <- c(15000,45000,80000,100000)
n_rep = n.rep = 10

#r2 mat represent the r2 matrix for the testing dataset
#column represent the ethnic groups
#row represent different p-value threshold

cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
setwd("/data/zhangh24/multi_ethnic/")


total <- 4*3
eth.vec <- rep("c",total)

l_vec <- rep(0,total)
m_vec <- rep(0,total)
ga_vc <- rep(0,total)
method_vec <- rep("c",total)
temp= 1
# for(i in 2:5){
#   for(i1 in 1:5){
#     for(l in 1:3){
#load the phenotype file
q = 1
result_list = list()
for(i in 2:5){
  for(l in 1:3){
    for(m in 1:4){
      for(i1 in 1:5){
        r2.test.rep <- rep(0,n.rep)
        y <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
        y <- y[,2+(1:n.rep)]
        n <- nrow(y)
        y_test_mat <- y[(100000+1):nrow(y),]
        #for(m in 1:4){
        
        
        
        for(i_rep in 1:n.rep){
          
          
          #for(k in 1:length(pthres)){
          
          #get the number of
          filename <- paste0(out.dir,eth[i],"/prs/prs_besteur_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".sscore")
          prs.temp <- as.data.frame(fread(filename))
          prs.score <- prs.temp[,4+q]
          y.test = y_test_mat[,i_rep]
          #y.vad = y_test_mat[(n.test+1):(nrow(y_test_mat)),i_rep]
          model1 <- lm(y.test~prs.score)
          r2.test.rep[i_rep] <- summary(model1)$r.square
          r2_function <- function(data, indices) {
            sample_data <- data[indices, ]
            model <- lm(y.vad ~ prs.vad, data = sample_data)
            return(summary(model)$r.square)
          }
          data <- data.frame(y.vad = y.test, prs.vad = prs.score)  # assuming y.vad and prs.vad are your data vectors
          # results <- boot(data = data, statistic = r2_function, R = 1000)
          # ci_result = boot.ci(results, type = "perc")
          # r2_low[i_rep] = ci_result$percent[4]
          # r2_high[i_rep] = ci_result$percent[5]
          
        }
        
        
        eursnp.result <- data.frame(eth.vec = rep(eth[i],n_rep),
                                    r2.vec = r2.test.rep,
                                    l_vec = rep(l,n_rep),
                                    m_vec = rep(m, n_rep),
                                    ga_vec= rep(i1, n_rep),
                                    rep_vec = c(1:n_rep),
                                    method_vec = rep(method[q],n_rep))
        result_list[[temp]] = eursnp.result
        temp = temp + 1
      }
    }
  }
}

eursnp.result = rbindlist(result_list)

#     }
#     }
#     
#     
#   }
# }


save(eursnp.result,file = paste0(out.dir,
                                 "eur.snp.reult.rep.rdata"))







#load(paste0(out.dir,"LD.clump.result.GA_",i1,".rdata"))


#save(r2.vec,file = paste0(out.dir,eth[i],"/r2_eursnp_",method[q],"_rho_",l,"_size_",m,"_GA_",i1))
#write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump.auc.csv")
