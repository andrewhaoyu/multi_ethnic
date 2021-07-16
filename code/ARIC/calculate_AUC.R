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
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")

qw = 2*3
eth_vec = rep("c",qw)
trait_vec = rep("c",qw)
r2_prs_vec = rep(0,qw)
rer2_prs_vec = rep(0,qw)
result.data.list = list()
reresult.data.list = list()

step = 1

for(i in 1:2){
  for(l in 1:3){
    
    temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
    data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
    out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
    
    
    y <- as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/pheno/pheno.txt")))
    colnames(y)[2] = "ID"
    # if(l==1){
    #   y$y = log(y$y)
    # }
    
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
    # sum.data = sum.data %>%
    #   mutate(BP=POS,SNP = SNP_ID,A1 = REF,
    #          P = PVAL)
    # idx.order = order(sum.data$CHR,sum.data$POS)
    # sum.data = sum.data[idx.order,] %>%
    #   mutate(index = c(1:nrow(sum.data)))
        clump.snp <- as.data.frame(fread(paste0(out.dir,"/LD_clump.clumped")))

    prs.clump <- left_join(clump.snp,sum.data,by="SNP")
    #prs.clump = prs.clump[order(prs.clump$index),]
    head(prs.clump)
    prs.filter = prs.clump %>% filter(CHR==1)
    tail(prs.filter)
    r2.vec.test.prs = rep(0,length(pthres))
    r2.vec.test.prs.pc = rep(0,length(pthres))
    r2.vec.vad.prs = rep(0,length(pthres))
    r2.vec.vad.prs.pc = rep(0,length(pthres))
    pthres_vec = rep(0,length(pthres))
    n.rep = 5
    
    
    
    r2.vec.test.prs.rep = matrix(0,length(pthres),n.rep)
    r2.vec.vad.prs.rep = matrix(0,length(pthres),n.rep)
    
    
    
    rer2.vec.test.prs.rep = matrix(0,length(pthres),n.rep)
    rer2.vec.vad.prs.rep = matrix(0,length(pthres),n.rep)
    
    
    
    
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
          #n.sub.fold = nrow(y)/n.rep
          start.end <- startend(nrow(y),n.rep,i_rep)
          vad.id = c(start.end[1]:start.end[2])
          test.id = setdiff(c(1:nrow(y)),vad.id)
          test.data <- y[test.id,]
          vad.data <- y[vad.id,]
          
          prs.test <- left_join(test.data,prs.data,by="ID")
          prs.vad <- left_join(vad.data,prs.data,by="ID")
          #model = lm(y~prs.score)
          
          model1.full <- lm(y~prs+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.test)
          model1.prs <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.test)
          #model1.null <- lm(y~age+sex,data=prs.test)
          #r2.test.rep[i_rep] <- summary(model1)$r.square
          model2.full <- lm(y~prs+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.vad)
          model2.prs <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.vad)
          #model2.null <- lm(y~age+sex,data=prs.vad)
          r2.vec.test.prs.rep[k,i_rep] = summary(model1.full)$r.square-summary(model1.prs)$r.square
          r2.vec.vad.prs.rep[k,i_rep] = summary(model2.full)$r.square-summary(model2.prs)$r.square
          
          #residual r2
          model1.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.test)
          model1.prs <- lm(model1.null$residual~prs,data=prs.test)
          #model1.null <- lm(y~age+sex,data=prs.test)
          #r2.test.rep[i_rep] <- summary(model1)$r.square
          model2.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.vad)
          model2.prs <- lm(model2.null$residual~prs,data=prs.vad)
          #model2.null <- lm(y~age+sex,data=prs.vad)
          rer2.vec.test.prs.rep[k,i_rep] = summary(model1.prs)$r.square
          rer2.vec.vad.prs.rep[k,i_rep] = summary(model2.prs)$r.square
          #pthres_vec[temp] = pthres[k]
          
          
        }
        
      }else{
        #pthres_vec[temp] = pthres[k]
        r2.vec.test.prs.rep[k,] = 0
        r2.vec.vad.prs.rep[k,] = 0
        rer2.vec.test.prs.rep[k,] = 0
        rer2.vec.vad.prs.rep[k,] = 0
        
      }
      
      
    }
    
    
    
    eth_vec[step] = eth[i]
    trait_vec[step] = trait[l]
    
    r2.vec.test.prs = rowMeans(r2.vec.test.prs.rep)
    r2.vec.vad.prs = rowMeans(r2.vec.vad.prs.rep)
    rer2.vec.test.prs = rowMeans(rer2.vec.test.prs.rep)
    rer2.vec.vad.prs = rowMeans(rer2.vec.vad.prs.rep)
    result.data <- data.frame(r2.vec.test.prs,r2.vec.vad.prs,
                              rer2.vec.test.prs,rer2.vec.vad.prs,
                              pthres,
                              eth = rep(eth[i],length(pthres)),
                              triat = rep(trait[l],length(pthres)))
    #find the best prs cutoff under each rep
    # r2.prs.temp = rep(0,n.rep) 
    # rer2.prs.temp = rep(0,n.rep)
    # for(i_rep in 1:n.rep){
    #   idx <- which.max(r2.vec.test.prs.rep[,i_rep]) 
    #   r2.prs.temp[i_rep] = r2.vec.vad.prs.rep[idx,i_rep]
    #   idx <- which.max(rer2.vec.test.prs.rep[,i_rep]) 
    #   rer2.prs.temp[i_rep] = rer2.vec.vad.prs.rep[idx,i_rep]
    # }
    # 
    # r2.prs = mean(r2.prs.temp)
    # 
    # rer2.prs = mean(rer2.prs.temp)
    
    idx <- which.max(result.data$r2.vec.test.prs)
    r2.prs = result.data$r2.vec.vad.prs[idx]
    idx <- which.max(result.data$rer2.vec.test.prs)
    rer2.prs = result.data$rer2.vec.vad.prs[idx]
    
    r2_prs_vec[step] = r2.prs
    rer2_prs_vec[step] = rer2.prs
    result.data.list[[step]] = result.data
    step = step+1
  }
  


r2.result = data.frame(eth = eth_vec,
                       trait = trait_vec,
                       r2_prs = r2_prs_vec,
                       rer2_prs = rer2_prs_vec
)
result.data = rbindlist(result.data.list)  
ARIC.result.CT = list(r2.result,result.data)
#save(ARIC.result.CT,file = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/ARIC.result.CT.rdata"))    
save(ARIC.result.CT,file = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/ARIC.result.CT.rep.rdata"))    

idx <- which(ARIC.result.CT[[2]]$eth =="EUR"&ARIC.result.CT[[2]]$triat =="urate")
ARIC.result.CT[[2]][idx,]
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
