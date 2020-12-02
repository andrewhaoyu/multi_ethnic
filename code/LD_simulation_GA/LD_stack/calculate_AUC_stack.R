#calculate AUC for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])
i1 = as.numeric(args[[5]])

library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
#n <- 120000
n.train.vec <- c(15000,45000,80000,100000)
#for(m in 1:1){
  #n.train <- n.train.vec[m]
  n.test <- 10000
  n.vad <- n.test
  n.rep = 10
  #r2 mat represent the r2 matrix for the testing dataset
  #column represent the ethnic groups
  #row represent different p-value threshold
  cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
  out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
  setwd("/data/zhangh24/multi_ethnic/")
  
  #load the phenotype file
  y <- as.data.frame(fread(paste0(cur.dir,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
  y <- y[,2+(1:n.rep),drop=F]
  n <- nrow(y)
  y_test_mat <- y[(100000+1):nrow(y),,drop=F]
  #y.test <- y[(n.train[i]+1):(n.train[i]+n.test[i])]
  # y.vad <- y[(n.train[i]+n.test[i]+1):n]
  #
  
  
  #for(i_rep in 1:n.rep){
    sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
    r2_vec = c(0.01,0.05,0.1,0.2,0.5)
    wc_base_vec = c(50,100,200,500)
    
    r2.vec.test <- rep(0,length(pthres)*length(r2_vec)*length(wc_base_vec))
    r2.vec.vad <- rep(0,length(pthres)*length(r2_vec)*length(wc_base_vec))
    pthres_vec <- rep(0,length(pthres)*length(r2_vec)*length(wc_base_vec))
    r2_ind_vec <- rep(0,length(pthres)*length(r2_vec)*length(wc_base_vec))
    wc_ind_vec <- rep(0,length(pthres)*length(r2_vec)*length(wc_base_vec))
    prs.mat <- matrix(0,n.test+n.vad,length(pthres)*length(r2_vec)*length(wc_base_vec))
    temp = 1
    for(r_ind in 1:length(r2_vec)){
      wc_vec = wc_base_vec/r2_vec[r_ind]
      for(w_ind in 1:length(wc_vec)){
        print(c(r_ind,w_ind))
        
    LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    clump.snp <- LD[,3,drop=F]  
    
    
    
    colnames(sum.data)[2] <- "SNP"
    #for(k in 1:length(pthres)){
    
    #for(k in 1:length(pthres)){
    
    prs.clump = left_join(clump.snp,sum.data,by="SNP")
    
    for(k in 1:length(pthres)){
      prs.all <- prs.clump %>% 
        filter(P<=pthres[k])
      if(nrow(prs.all)>0){
            
            
            filename <- paste0(out.dir,eth[i],"/prs/prs_pvalue_",k,"_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".profile")
            
            prs.temp <- fread(filename)
            prs.score <- prs.temp$V1
            #prs.score <- prs.temp$SCORE*2*length(idx)+prs.score
        #prs.score.mat[,k] = prs.score
        prs.test <- prs.score[(1):(n.test)]
        prs.vad <- prs.score[(n.test+1):(n.test+n.vad)]
        #model = lm(y~prs.score)
        y.test = y_test_mat[1:n.test,i_rep]
        y.vad = y_test_mat[(n.test+1):(nrow(y_test_mat)),i_rep]
        model1 <- lm(y.test~prs.test)
        #r2.test.rep[i_rep] <- summary(model1)$r.square
        model2 <- lm(y.vad~prs.vad)
        r2.vec.test[temp] = summary(model1)$r.square
        r2.vec.vad[temp] = summary(model2)$r.square
        pthres_vec[temp] = pthres[k]
        r2_ind_vec[temp] = r_ind
        wc_ind_vec[temp] = w_ind
        prs.mat[,temp] = prs.score
        temp = temp+1
      }else{
        r2.vec.test[temp] = 0
        r2.vec.vad[temp] = 0  
        pthres_vec[temp] = pthres[k]
        r2_ind_vec[temp] = r_ind
        wc_ind_vec[temp] = w_ind
        prs.mat[,temp] = 0
        temp = temp+1
      }
      
    }
  }
    
     }
  
    # }
    
    
    #}
    prs.sum = colSums(prs.mat)
    idx <- which(prs.sum!=0)
    #drop the prs with all 0
    prs.mat <- prs.mat[,idx]
    library(SuperLearner)
    library(ranger)
    x.test = as.data.frame(prs.mat[1:n.test,])
    x.vad= as.data.frame(prs.mat[(1+n.test):(n.test+n.vad),])
    SL.libray <- c(
                  #"SL.xgboost"
                   #"SL.randomForest"
                  "SL.glmnet",
                  "SL.ridge",
                  #"SL.bayesglm"
                  #"SL.stepAIC"
                   "SL.nnet"
                   #"SL.ksvm",
                  #"SL.bartMachine", 
                 #"SL.kernelKnn",
                #"SL.rpartPrune", 
                 #"SL.lm"
                #"SL.mean"
      )
    sl = SuperLearner(Y = y.test, X = x.test, family = gaussian(),
                            # For a real analysis we would use V = 10.
                           # V = 3,
                            SL.library = SL.libray)
    sl
    y.pred <- predict(sl, x.vad, onlySL = TRUE)
    #names(r2.vec.test) <- names(r2.vec.vad) <- pthres
    
    #evaluate the best prs performance on the validation
    model <- lm(y.vad~y.pred[[1]])
    r2.stack <- summary(model)$r.square
    result.data <- data.frame(r2.vec.test,r2.vec.vad,
                              pthres_vec,r2_ind_vec,
                              wc_ind_vec)
    #standard C+T
    result.data.CT = result.data %>% 
      filter(r2_ind_vec==3&wc_ind_vec==1)
    idx <- which.max(result.data.CT$r2.vec.test)
    r2.max.ct <- result.data.CT$r2.vec.vad[idx]
    idx <- which.max(r2.vec.test)
    r2.max <- r2.vec.vad[idx]
    r2.list <- list(r2.stack,
                    r2.max,
                    r2.max.ct)
    save(r2.list,file = paste0(out.dir,eth[i],"/r2.list_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
  #}
#}
#write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump.auc.csv")
