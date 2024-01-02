#calculate the AUC for EUR SNPs with EUR coefficients
#with target ethnic coefficients
#with EB coefficients
#calculate AUC for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
#q for three different methods

#i_rep = as.numeric(args[[5]])
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
i1 = as.numeric(args[[2]])
l = as.numeric(args[[3]])
m = as.numeric(args[[4]])
library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
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
# result_list = list()
# temp_list = 1
# for(i in 2:5){
#   for(i1 in 1:5){
#     for(l in 1:3){
#       for(m in 1:4){
        
        total <- 1
        eth.vec <- rep("c",total)
        r2.result <- ci_low <- ci_high <- rep(0,total)
        l_vec <- rep(0,total)
        m_vec <- rep(0,total)
        ga_vc <- rep(0,total)
        method_vec <- rep("PRS-CSx (five ancestries)",total)
        temp= 1
        #phi = c("1e-02","1e-04","1e-06")
        phi = c("1e+00","1e-02","1e-04","1e-06")
        # for(i in 2:5){
        #   for(i1 in 1:5){
        #     for(l in 1:3){
        #load the phenotype file
        y <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
        y <- y[,2+(1:n.rep)]
        n <- nrow(y)
        y_test_mat <- y[(100000+1):nrow(y),]
        #for(m in 1:4){
        #for(q in 1:3){
        
        print(m)
        
        r2.vad.rep = r2_low = r2_high <- rep(0,n.rep)
        
        for(i_rep in 1:n.rep){
          r2.test.rep <- rep(0,3)
          coef.mat = matrix(NA,4,length(eth))
          #find the best phi
          for(v in 1:4){
            prs_list = list()
            #load target prs
            for(i_eth in 1:5){
              filename <- paste0(out.dir.sum,eth[i],"/prscsx_five/prs_csx_",eth[i_eth],"_rho_",l,"_size_",m,"_GA_",i1,"_phi",phi[v],".sscore")  
              prs_file = as.data.frame(fread(filename))
              prs_list[[i_eth]] = prs_file[,4+i_rep]
            }
            prs_file = bind_cols(prs_list)
            prs_tun = as.matrix(prs_file[1:10000,])
            y.test = y_test_mat[1:10000,i_rep]
            
            model1 <- lm(y.test~prs_tun)
            coef.mat[v,] = coefficients(model1)[2:(ncol(prs_tun)+1)]
            r2.test.rep[v] <- summary(model1)$r.square
            
          }
          idx <- which.max(r2.test.rep)
          #load the max prs
          prs_list = list()
          #load target prs
          for(i_eth in 1:5){
            filename <- paste0(out.dir.sum,eth[i],"/prscsx_five/prs_csx_",eth[i_eth],"_rho_",l,"_size_",m,"_GA_",i1,"_phi",phi[idx],".sscore")  
            prs_file = as.data.frame(fread(filename))
            prs_list[[i_eth]] = prs_file[,4+i_rep]
          }
          prs_file = bind_cols(prs_list)
          prs_vad = as.matrix(prs_file[10001:20000,])
          prs = prs_vad%*%coef.mat[idx,]
          y.vad = y_test_mat[10001:20000,i_rep]
          model = lm(y.vad~prs)
          r2.vad.rep[i_rep] <- summary(model)$r.square
        }
        
        
        eth.vec <- rep(eth[i],n.rep)
        l_vec <- rep(l,n.rep)
        m_vec <- rep(m,n.rep)
        ga_vc <- rep(i1,n.rep)
        #method_vec[temp] <- method[q]
        
        
        #}
        # }
        #     }
        #     
        #     
        #   }
        # }
        
        prscsx.result <- data.frame(eth.vec,
                                    r2.vec = r2.vad.rep,
                                    l_vec,
                                    m_vec,
                                    ga_vec=ga_vc,
                                    rep_vec = c(1:n.rep),
                                    method_vec = method_vec)
#         result_list[[temp_list]] = prscsx.result
#         temp_list = temp_list + 1
#       }
#     }
#   }
# }

save(prscsx.result,file = paste0(out.dir,
                                 "prscsx_five.result_rep_ga_",i1,"_eth_",i,"_rho_",l,
                                 "_size_",m,".rdata"))




