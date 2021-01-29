startend <- function(num,size,ind){
  split.all <- split(1:num,cut(1:num,size))
  temp <- split.all[[ind]]
  start <- temp[1]
  end <- temp[length(temp)]
  return(c(start,end))
}
#library(devtools)
#install.packages("devtools")
library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")

#one ethnic group
#three trait
qw = 1*3
eth_vec = rep("c",qw)
trait_vec = rep("c",qw)
rer2_prs_vec = rep(0,qw)
rer2_prs_sl<- rep(0,qw)
result.data.list = list()
reresult.data.list = list()

#step for loop in trait and eth
step = 1

i = 2
  for(l in 1:3){
    
    temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
    data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
    out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
    
    y <- as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/pheno/pheno.txt")))
    colnames(y)[2] = "ID"
    
    
    n.col = length(pthres)^2*length(r2_vec)*length(wc_base_vec)
    
    r2.vec.test.prs = rep(0,n.col)
    r2.vec.test.prs.pc = rep(0,n.col)
    r2.vec.vad.prs = rep(0,n.col)
    r2.vec.vad.prs.pc = rep(0,n.col)
    prs.mat = matrix(0,)
    n.rep = 5
    
    genotype.fam <- as.data.frame(fread(paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc1.fam")))
    
    n.sub = nrow(genotype.fam)
    r2.vec.test.prs.rep = matrix(0,n.col,n.rep)
    r2.vec.vad.prs.rep = matrix(0,n.col,n.rep)
    rer2.vec.test.prs.rep = matrix(0,n.col,n.rep)
    rer2.vec.vad.prs.rep = matrix(0,n.col,n.rep)
    prs.mat = data.frame(ID = genotype.fam$V2,matrix(0,n.sub,n.col))
    #files = dir(path = temp.dir,pattern=paste0(".profile"),full.names = T)
    #temp for loop in rbind,w_ind,k1,k2
    #evaluate 2DLD-method result
    temp =1
    for(r_ind in 1:length(r2_vec)){
      wc_vec = wc_base_vec/r2_vec[r_ind]
      for(w_ind in 1:length(wc_vec)){
        print(c(r_ind,w_ind))
        
    
    for(k1 in 1:length(pthres)){
      for(k2 in 1:length(pthres)){
        
      
        filename <- paste0(temp.dir,"/2Dprs_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,"_",k2,".profile")
        
        prs.temp <- fread(filename)
        prs.score <- prs.temp$V1
        
        
        prs.data = data.frame(ID = genotype.fam$V2,prs= prs.score,stringsAsFactors = F)
        prs.mat[,temp+1] = prs.data$prs
        
        #prs.score <- prs.temp$SCORE*2*length(idx)+prs.score
        #prs.score.mat[,k] = prs.score
        for(i_rep in 1:n.rep){
          #n.sub.fold = nrow(y)/n.rep
          start.end <- startend(nrow(y),n.rep,i_rep)
          vad.id = c(start.end[1]:start.end[2])
          test.id = setdiff(c(1:nrow(y)),vad.id)
          test.data <- y[test.id,]
          vad.data <- y[vad.id,]
          
          prs.test <- left_join(test.data,prs.data,by="ID")
          prs.vad <- left_join(vad.data,prs.data,by="ID")
          #model = lm(y~prs.score)
          
          # model1.full <- lm(y~prs+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.test)
          # model1.prs <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.test)
          # #model1.null <- lm(y~age+sex,data=prs.test)
          # #r2.test.rep[i_rep] <- summary(model1)$r.square
          # model2.full <- lm(y~prs+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.vad)
          # model2.prs <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.vad)
          # #model2.null <- lm(y~age+sex,data=prs.vad)
          # r2.vec.test.prs.rep[temp,i_rep] = summary(model1.full)$r.square-summary(model1.prs)$r.square
          # r2.vec.vad.prs.rep[temp,i_rep] = summary(model2.full)$r.square-summary(model2.prs)$r.square
          
          #residual r2
          model1.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.test)
          model1.prs <- lm(model1.null$residual~prs,data=prs.test)
          #model1.null <- lm(y~age+sex,data=prs.test)
          #r2.test.rep[i_rep] <- summary(model1)$r.square
          model2.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.vad)
          model2.prs <- lm(model2.null$residual~prs,data=prs.vad)
          #model2.null <- lm(y~age+sex,data=prs.vad)
          rer2.vec.test.prs.rep[temp,i_rep] = summary(model1.prs)$r.square
          rer2.vec.vad.prs.rep[temp,i_rep] = summary(model2.prs)$r.square
          #pthres_vec[temp] = pthres[k]
          
          
        }
        temp = temp+1
      }
      } 
     
        
        
        
         
    }
    }
    eth_vec[step] = eth[i]
    trait_vec[step] = trait[l]
    
    # r2.vec.test.prs = rowMeans(r2.vec.test.prs.rep)
    # r2.vec.vad.prs = rowMeans(r2.vec.vad.prs.rep)
    rer2.vec.test.prs = rowMeans(rer2.vec.test.prs.rep)
    rer2.vec.vad.prs = rowMeans(rer2.vec.vad.prs.rep)
    result.data <- data.frame(rer2.vec.test.prs,rer2.vec.vad.prs,
                              expand.grid(pthres,pthres)[,c(2,1)],
                              eth = rep(eth[i],length(pthres)^2),
                              triat = rep(trait[l],length(pthres)^2))
    
    
    idx <- which.max(result.data$r2.vec.test.prs)
    r2.prs = result.data$r2.vec.vad.prs[idx]
    idx <- which.max(result.data$rer2.vec.test.prs)
    rer2.prs = result.data$rer2.vec.vad.prs[idx]
    
    r2_prs_vec[step] = r2.prs
    rer2_prs_vec[step] = rer2.prs
    result.data.list[[step]] = result.data
    
    #evaluate 2DLD-super-learning method
    
    #r2_prs_sl <- rep(0,qw)
    
    #r2.prs.sl.rep <- rep(0,n.rep)
    rer2.prs.sl.rep <- rep(0,n.rep)
    SL.libray <- c(
      "SL.glmnet",
      "SL.ridge"
      #"SL.nnet"
    )
    library(caret)
    library(SuperLearner)
    library(ranger)
    library(glmnet)
    for(i_rep in 1:n.rep){
      #n.sub.fold = nrow(y)/n.rep
      start.end <- startend(nrow(y),n.rep,i_rep)
      vad.id = c(start.end[1]:start.end[2])
      test.id = setdiff(c(1:nrow(y)),vad.id)
      test.data <- y[test.id,]
      vad.data <- y[vad.id,]
      
      prs.test <- left_join(test.data,prs.mat,by="ID")
      prs.vad <- left_join(vad.data,prs.mat,by="ID")
      #model = lm(y~prs.score)
      
      #residual r2
      model1.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.test)
      y.test <- model1.null$residual
      x.test = prs.test[,16:ncol(prs.test)]
      sl = SuperLearner(Y = y.test, X = x.test, family = gaussian(),
                        # For a real analysis we would use V = 10.
                        # V = 3,
                        SL.library = SL.libray)
      x.vad = prs.vad[,16:ncol(prs.vad)]
      y.pred <- predict(sl,x.vad,onlySL = TRUE)
      
      model2.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.vad)
      y.vad <- model2.null$residual
      
      model2.prs <- lm(y.vad~y.pred[[1]])
      #model2.null <- lm(y~age+sex,data=prs.vad)
      rer2.prs.sl.rep[i_rep] = summary(model2.prs)$r.square
      #pthres_vec[temp] = pthres[k]
      
      
    }
   
    rer2_prs_sl = mean(rer2_prs_sl)
    step = step+1
    
  }  

r2.result = data.frame(eth = eth_vec,
                       trait = trait_vec,
                       rer2_prs = rer2_prs_vec,
                       rer2_prs_sl = rer2_prs_sl
)
result.data = rbindlist(result.data.list)  
ARIC.result.2DLD = list(r2.result,result.data)
#save(ARIC.result.CT,file = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/ARIC.result.CT.rdata"))    
save(ARIC.result.2DLD,file = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/ARIC.result.2DLD.rdata"))    
#}
#     prs.sum = colSums(prs.mat)
#     idx <- which(prs.sum!=0)
#     #drop the prs with all 0
#     prs.mat <- prs.mat[,idx]
#     library(SuperLearner)
#     library(ranger)
#     x.test = as.data.frame(prs.mat[1:n.test,])
#     x.vad= as.data.frame(prs.mat[(1+n.test):(n.test+n.vad),])
#     SL.libray <- c(
#       #"SL.xgboost"
#       #"SL.randomForest"
#       "SL.glmnet",
#       "SL.ridge",
#       #"SL.bayesglm"
#       #"SL.stepAIC"
#       "SL.nnet"
#       #"SL.ksvm",
#       #"SL.bartMachine", 
#       #"SL.kernelKnn",
#       #"SL.rpartPrune", 
#       #"SL.lm"
#       #"SL.mean"
#     )
#     sl = SuperLearner(Y = y.test, X = x.test, family = gaussian(),
#                       # For a real analysis we would use V = 10.
#                       # V = 3,
#                       SL.library = SL.libray)
#     sl
#     y.pred <- predict(sl, x.vad, onlySL = TRUE)
#     #names(r2.vec.test) <- names(r2.vec.vad) <- pthres
#     
#     #evaluate the best prs performance on the validation
#     model <- lm(y.vad~y.pred[[1]])
#     r2.stack <- summary(model)$r.square
#     result.data <- data.frame(r2.vec.test,r2.vec.vad,
#                               pthres_vec,r2_ind_vec,
#                               wc_ind_vec)
#     #standard C+T
#     result.data.CT = result.data %>% 
#       filter(r2_ind_vec==3&wc_ind_vec==1)
#     idx <- which.max(result.data.CT$r2.vec.test)
#     r2.max.ct <- result.data.CT$r2.vec.vad[idx]
#     idx <- which.max(r2.vec.test)
#     r2.max <- r2.vec.vad[idx]
#     r2.list <- list(r2.stack,
#                     r2.max,
#                     r2.max.ct)
#     save(r2.list,file = paste0(out.dir,eth[i],"/r2.list_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
#     
#   }
# }

#}
#}
#write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump.auc.csv")
