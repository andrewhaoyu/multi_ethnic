#merge the prs by chromosome to one for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
args = commandArgs(trailingOnly = T)
i = 2
l = as.numeric(args[[1]])
m = as.numeric(args[[2]])
i_rep = as.numeric(args[[3]])
#i_rep = 2
i1 = 1

library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
#n <- 120000

#for(m in 1:1){

n.test <- 10000
n.vad <- n.test
n.rep = 10
#r2 mat represent the r2 matrix for the testing dataset
#column represent the ethnic groups
#row represent different p-value threshold
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
setwd("/data/zhangh24/multi_ethnic/")
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
y <- as.data.frame(fread(paste0(cur.dir,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
y <- y[,2+(1:n.rep),drop=F]
n <- nrow(y)
y_test_mat <- y[(100000+1):nrow(y),,drop=F]

#for(i_rep in 1:n.rep){
summary.eur <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))
colnames(summary.eur)[9] = "peur"
colnames(summary.eur)[7] = "beta_eur"
summary.eur.select = summary.eur %>%
  select(SNP,beta_eur,peur)
sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))
colnames(sum.data)[2] <- "SNP"
#combine the target level summary stat with EUR
summary.com <- left_join(sum.data,summary.eur.select,by="SNP")

r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)

r2.vec.test <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
r2.vec.vad <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
pthres_vec1 <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
pthres_vec2 <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
r2_ind_vec <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
wc_ind_vec <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
prs.mat <- matrix(0,n.test+n.vad,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
temp = 1

for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    
    LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    clump.snp <- LD
    
    #    colnames(sum.data)[2] <- "SNP"
    
    #for(k in 1:length(pthres)){
    
    prs.clump = left_join(clump.snp,summary.com,by="SNP")
    
    for(k1 in 1:length(pthres)){
      for(k2 in 1:length(pthres)){
        prs.all <- prs.clump %>%
          filter(peur<=pthres[k1]|
                   P<=pthres[k2])
        if(nrow(prs.all)>0){
          filename <- paste0(out.dir,eth[i],"/prs/prs_eb_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,".p_value_",k2,".profile")
          
          prs.temp <- fread(filename)  
          prs.score <- prs.temp$SCORE
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
          pthres_vec1[temp] = pthres[k1]
          pthres_vec2[temp] = pthres[k2]
          r2_ind_vec[temp] = r_ind
          wc_ind_vec[temp] = w_ind
          prs.mat[,temp] = prs.score
          temp = temp+1
        }else{
          r2.vec.test[temp] = 0
          r2.vec.vad[temp] = 0
          pthres_vec1[temp] = pthres[k1]
          pthres_vec2[temp] = pthres[k2]
          r2_ind_vec[temp] = r_ind
          wc_ind_vec[temp] = w_ind
          prs.mat[,temp] = 0
          temp = temp+1
        }
        
      }
    }
    
  }
}
result.data <- data.frame(r2.vec.test,r2.vec.vad,
                          pthres_vec1,pthres_vec2,
                          r2_ind_vec,
                          wc_ind_vec)
prs.sum = colSums(prs.mat)
idx <- which(prs.sum!=0)
#drop the prs with all 0
prs.mat <- prs.mat[,idx]
#drop the columns with perfect correlation
prs.mat = as.data.frame(prs.mat)
mtx = cor(prs.mat[1:n.test,])
library(caret)
drop = findCorrelation(mtx,cutoff=0.98)
drop = names(prs.mat)[drop]
prs.mat.new = prs.mat %>% 
  select(-all_of(drop))
library(SuperLearner)
library(ranger)
x.test = as.data.frame(prs.mat.new[1:n.test,])
x.vad= as.data.frame(prs.mat.new[(1+n.test):(n.test+n.vad),])

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
y.pred.test = predict(sl,x.test,onlySL = TRUE)
model <- lm(y.test~y.pred.test[[1]])
r2.stack.test <- summary(model)$r.square

y.pred <- predict(sl, x.vad, onlySL = TRUE)
#names(r2.vec.test) <- names(r2.vec.vad) <- pthres

#evaluate the best prs performance on the validation
model <- lm(y.vad~y.pred[[1]])
r2.stack <- summary(model)$r.square
result.data <- data.frame(r2.vec.test,r2.vec.vad,
                          pthres_vec1,pthres_vec2,
                          r2_ind_vec,
                          wc_ind_vec)
#standard C+T
idx <- which.max(r2.vec.test)
r2.max <- r2.vec.vad[idx]
r2.list = c(r2.stack.test,r2.stack)
save(r2.list,file = paste0(out.dir,eth[i],"/r2_overfit_rho_eb_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
#write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump.auc.csv")
