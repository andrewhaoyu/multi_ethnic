

library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")

qw = 3*3
method_vec = rep("c",qw)
eth_vec = rep("c",qw)
trait_vec = rep("c",qw)
r2_prs_vec = rep(0,qw)
rer2_prs_pc_vec = rep(0,qw)
result.data.list = list()
method_op <- c("Best EUR SNP (C+T)",
               "Best EUR SNP + target coefficients (C+T)",
               "Best EUR SNP + EB coefficients (C+T)")
               
step = 1
i= 2
#for(i in 1:2){
  for(l in 1:3){
    
    temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
    data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
    out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
    
    y <- as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/pheno/pheno.txt")))
    colnames(y)[2] = "ID"
    sum.data = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/sumdata/training-GWAS-formatted.txt")))
    
    sum.data = sum.data %>% 
      mutate(BP=POS,SNP = SNP_ID,A1 = REF,
             P = PVAL) 
        
        filename <- paste0(out.dir,"prs_best_eur_prs.profile")
        
        prs.temp <- as.data.frame(fread(filename,header=T))
        for(i_method in 1:3){
          prs.score <-   prs.temp[,i_method]
          genotype.fam <- as.data.frame(fread(paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc1.fam")))
          prs.data = data.frame(ID = genotype.fam$V2,prs= prs.score,stringsAsFactors = F)
          #prs.score <- prs.temp$SCORE*2*length(idx)+prs.score
          #prs.score.mat[,k] = prs.score
          prs.test <- left_join(y,prs.data,by="ID")
          
          model1.full <- lm(y~prs+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.test)
          model1.prs <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.test)
          #model1.null <- lm(y~age+sex,data=prs.test)
          #r2.test.rep[i_rep] <- summary(model1)$r.square
          
          
          r2_prs_vec[step] = summary(model1.full)$r.square-summary(model1.prs)$r.square
          model1.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.test)
          model1.prs <- lm(model1.null$residual~prs,data=prs.test)
          rer2_prs_pc_vec[step] = summary(model1.prs)$r.square
          method_vec[step] = method_op[i_method]
          eth_vec[step] = eth[i]
          trait_vec[step] = trait[l]
          
          step = step+ 1
        }
        
        
  }
      
    
  
#}

r2.result = data.frame(eth = eth_vec,
                       trait = trait_vec,
                       r2_prs = r2_prs_vec,
                       rer2_prs_pc = rer2_prs_pc_vec,
                       method_vec)

ARIC.result.bestEUR = r2.result
save(ARIC.result.bestEUR,file = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/ARIC.result.bestEUR.rdata"))    
write.csv(ARIC.result.bestEUR,file = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/ARIC.result.bestEUR.csv"))

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