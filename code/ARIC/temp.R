

library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")

qw = 2*3
eth_vec = rep("c",qw)
trait_vec = rep("c",qw)
r2_prs_vec = rep(0,qw)
r2_prs_pc_vec = rep(0,qw)
result.data.list = list()

step = 1

for(i in 1:2){
  for(l in 1:3){
    
    temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
    data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
    out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
    
    y <- as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/pheno/pheno.txt")))
    colnames(y)[2] = "ID"
    # test.id = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/pheno/IDs.tuning.txt"),header=F))
    # colnames(test.id)[2] = "ID"
    # vad.id = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/pheno/IDs.validation.txt"),header=F))
    # colnames(vad.id)[2] = "ID"
    # test.id = sample(c(1:nrow(y)),nrow(y)/2,replace = F)
    # test.data = left_join(test.id,y,by="ID")
    # vad.data = left_join(vad.id,y,by="ID")
    sum.data = as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/sumdata/training-GWAS-formatted.txt")))
    
    sum.data = sum.data %>% 
      mutate(BP=POS,SNP = SNP_ID,A1 = REF,
             P = PVAL) 
    clump.snp <- as.data.frame(fread(paste0(out.dir,"/LD_clump.clumped")))
    
    prs.clump <- left_join(clump.snp,sum.data,by="SNP")
    r2.vec.test.prs = rep(0,length(pthres))
    r2.vec.test.prs.pc = rep(0,length(pthres))
    r2.vec.vad.prs = rep(0,length(pthres))
    r2.vec.vad.prs.pc = rep(0,length(pthres))
    pthres_vec = rep(0,length(pthres))
    n.rep = 100
    
    
    
    r2.vec.test.prs.rep = rep(0,n.rep)
    r2.vec.test.prs.pc.rep = rep(0,n.rep)
    r2.vec.vad.prs.rep = rep(0,n.rep)
    r2.vec.vad.prs.pc.rep = rep(0,n.rep)
    
    temp =1
    for(k in 1:length(pthres)){
      prs.all <- prs.clump %>% 
        filter(P<=pthres[k])
      if(nrow(prs.all)>0){
        
        
        filename <- paste0(out.dir,"/prs_pvalue_",k,".profile")
        
        prs.temp <- fread(filename)
        prs.score <- prs.temp$V1
        genotype.fam <- as.data.frame(fread(paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc1.fam")))
        
        prs.data = data.frame(ID = genotype.fam$V2,prs= prs.score,stringsAsFactors = F)
        #prs.score <- prs.temp$SCORE*2*length(idx)+prs.score
        #prs.score.mat[,k] = prs.score
        for(i_rep in 1:n.rep){
          test.id = sample(c(1:nrow(y)),nrow(y)/2,replace = F)
          vad.id = setdiff(c(1:nrow(y)),test.id)
          test.data <- y[test.id,]
          vad.data <- y[vad.id,]
          
          prs.test <- left_join(test.data,prs.data,by="ID")
          prs.vad <- left_join(vad.data,prs.data,by="ID")
          #model = lm(y~prs.score)
          
          model1.full <- lm(y~prs+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.test)
          model1.prs <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.test)
          model1.null <- lm(y~age+sex,data=prs.test)
          #r2.test.rep[i_rep] <- summary(model1)$r.square
          model2.full <- lm(y~prs+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.vad)
          model2.prs <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex+age+sex,data=prs.vad)
          model2.null <- lm(y~age+sex,data=prs.vad)
          r2.vec.test.prs.rep[i_rep] = summary(model1.full)$r.square-summary(model1.prs)$r.square
          r2.vec.vad.prs.rep[i_rep] = summary(model2.full)$r.square-summary(model2.prs)$r.square
          r2.vec.test.prs.pc.rep[i_rep] = summary(model1.full)$r.square-summary(model1.null)$r.square
          r2.vec.vad.prs.pc.rep[i_rep] = summary(model2.full)$r.square-summary(model2.null)$r.square
          
          #coefficients(model1.full)['scale(prs)']^2/var(prs.test$y)
          
          idx <- which(is.na(prs.test$y))
          length(idx)            
          
          # r2.vec.test.prs[temp] = summary(model1.prs)$r.square
          # r2.vec.vad.prs[temp] = summary(model2.prs)$r.square
          # r2.vec.test.prs.pc[temp] = summary(model1.full)$r.square
          # r2.vec.vad.prs.pc[temp] = summary(model2.full)$r.square
          # 
        }
        pthres_vec[temp] = pthres[k]
        r2.vec.test.prs[temp] = mean(r2.vec.test.prs.rep)
        r2.vec.vad.prs[temp] = mean(r2.vec.vad.prs.rep)
        r2.vec.test.prs.pc[temp] = mean(r2.vec.test.prs.pc.rep)
        r2.vec.vad.prs.pc[temp] = mean(r2.vec.vad.prs.pc.rep)
        temp = temp+1 
      }
      else{
        pthres_vec[temp] = pthres[k]
        r2.vec.test.prs[temp] = 0
        r2.vec.vad.prs[temp] = 0
        r2.vec.test.prs.pc[temp] = 0
        r2.vec.vad.prs.pc[temp] = 0
        temp = temp+1 
      }
      
      
      
    }
    
    
    
    eth_vec[step] = eth[i]
    trait_vec[step] = trait[l]
    result.data <- data.frame(r2.vec.test.prs,r2.vec.vad.prs,
                              r2.vec.test.prs.pc,r2.vec.vad.prs.pc,
                              pthres_vec,
                              eth = rep(eth[i],length(pthres)),
                              triat = rep(trait[l],length(pthres)))
    
    idx <- which.max(r2.vec.test.prs)
    r2.prs = r2.vec.vad.prs[idx]
    idx <- which.max(r2.vec.test.prs.pc)
    r2.prs.pc = r2.vec.vad.prs.pc[idx]
    
    r2_prs_vec[step] = r2.prs
    r2_prs_pc_vec[step] = r2.prs.pc
    result.data.list[[step]] = result.data
    step = step+1
  }
  
}

r2.result = data.frame(eth = eth_vec,
                       trait = trait_vec,
                       r2_prs = r2_prs_vec,
                       r2_prs_pc = r2_prs_pc_vec)
result.data = rbindlist(result.data.list)  
ARIC.result.CT = list(r2.result,result.data)
#save(ARIC.result.CT,file = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/ARIC.result.CT.rdata"))    
save(ARIC.result.CT,file = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/ARIC.result.CT.rep.rdata"))    
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
