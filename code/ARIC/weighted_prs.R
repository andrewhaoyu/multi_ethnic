#weighted prs of AFR and best EUR

library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")

qw = 3
method_vec = rep("c",qw)
eth_vec = rep("c",qw)
trait_vec = rep("c",qw)
r2_prs_vec = rep(0,qw)
rer2_prs_pc_vec = rep(0,qw)
result.data.list = list()
# method_op <- c("Best EUR SNP (C+T)",
#                "Best EUR SNP + target coefficients (C+T)",
#                "Best EUR SNP + EB coefficients (C+T)")

step = 1
i= 2
i_method = 1
#for(i in 1:2){
for(l in 1:3){
  
  temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
  data.dir = "/dcl01/chatterj/data/jin/prs/realdata/ARIC/"
  out.dir = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/",trait[l],"/",eth[i],"/")
  
  y <- as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/pheno/pheno.txt")))
  colnames(y)[2] = "ID"
  
  filename <- paste0(out.dir,"prs_best_eur_prs.profile")
  
  prs.temp <- as.data.frame(fread(filename,header=T))
  prs.eur = prs.temp[,i_method]
  
  #load best target ethnic group prs
  load("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/ARIC.result.CT.rep.rdata")
  ARIC.result.CT = ARIC.result.CT[[2]]
  colnames(ARIC.result.CT)[7] = "trait_vec"
  ARIC.result.CT.temp =  ARIC.result.CT %>%
    filter(trait_vec ==trait[l]&
             eth == "AFR") %>% 
    filter(rer2.vec.test.prs==max(rer2.vec.test.prs))
  k = which(pthres==ARIC.result.CT.temp$pthres)
  filename <- paste0(out.dir,"/prs_pvalue_",k,".profile")
  prs.temp <- fread(filename)
  prs.tar <- prs.temp$V1
  
  
    prs.score <-   cbind(prs.eur,prs.tar)
    genotype.fam <- as.data.frame(fread(paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc1.fam")))
    prs.data = data.frame(ID = genotype.fam$V2, prs.score,stringsAsFactors = F)
    #prs.score <- prs.temp$SCORE*2*length(idx)+prs.score
    #prs.score.mat[,k] = prs.score
    
    prs.all <- left_join(y,prs.data,by="ID")
    
    #
    startend <- function(num,size,ind){
      split.all <- split(1:num,cut(1:num,size))
      temp <- split.all[[ind]]
      start <- temp[1]
      end <- temp[length(temp)]
      return(c(start,end))
    }
    n.rep = 5
    rer2_prs_pc_rep =rep(0,n.rep)
    #split the data into 5 fold
    for(i_rep in 1:n.rep){
      start.end <- startend(nrow(prs.all),n.rep,i_rep)
      vad.id = c(start.end[1]:start.end[2])
      test.id = setdiff(c(1:nrow(prs.all)),vad.id)
      
      prs.test = prs.all[test.id,]
      prs.vad = prs.all[vad.id,]
      
      model1.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.test)
      model1.prs <- lm(model1.null$residual~prs.tar+prs.eur,data=prs.test)
      
      beta = coefficients(model1.prs)
      
      prs.vad.score <- cbind(prs.vad$prs.tar,prs.vad$prs.eur)%*%beta[2:3]
      model2.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=prs.vad)
      model2.prs <- lm(model2.null$residual~prs.vad.score)
      rer2_prs_pc_rep[i_rep] = summary(model2.prs)$r.square
    }
    rer2_prs_pc_vec[step] = mean(rer2_prs_pc_rep)
    method_vec[step] = "Weighted PRS"
    eth_vec[step] = eth[i]
    trait_vec[step] = trait[l]
    
    step = step+ 1
  }
  
  




#}

weightedprs.result = data.frame(eth = eth_vec,
                       trait = trait_vec,
                       r2_prs = r2_prs_vec,
                       rer2_prs = rer2_prs_pc_vec,
                       method_vec)
save(weightedprs.result,file = "/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/ARIC.weightedprs.result.rdata")

# ARIC.result.bestEUR = r2.result
# save(ARIC.result.bestEUR,file = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/ARIC.result.bestEUR.rdata"))    
# write.csv(ARIC.result.bestEUR,file = paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/multi_ethnic/result/ARIC/ARIC.result.bestEUR.csv"))
