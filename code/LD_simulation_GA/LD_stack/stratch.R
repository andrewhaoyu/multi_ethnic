#merge the prs by chromosome to one for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])
#i_rep = 2
i1 = as.numeric(args[[5]])

library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
#n <- 120000

#for(m in 1:1){

n.test <- 10000
n.vad <- n.test
n.rep = 3
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

# r2_vec = c(0.01,0.05,0.1,0.2,0.5)
# wc_base_vec = c(50,100,200,500)
r2_vec = c(0.01,0.05,0.1,0.2,0.5)
wc_base_vec = c(50)

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
    clump.snp <- LD[,3,drop=F]  
    
    colnames(sum.data)[2] <- "SNP"
    
    #for(k in 1:length(pthres)){
    
    prs.clump = left_join(clump.snp,summary.com,by="SNP")
    
    for(k1 in 1:length(pthres)){
      for(k2 in 1:length(pthres)){
        prs.all <- prs.clump %>% 
          filter(peur<=pthres[k1]|
                   P<=pthres[k2])
        if(nrow(prs.all)>0){
          filename <- paste0(out.dir,eth[i],"/prs/prs_eb_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,"p_value,",k2,".profile")
          
          prs.temp <- fread(filename)  
          prs.score <- prs.temp$V1
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

pthres_new <- c(5E-08,5E-07,5E-06,5E-05,1E-04,1E-03,1E-02,0.5)
result.data = result.data %>% 
  filter(pthres_vec1%in%pthres_new&
           pthres_vec2%in%pthres_new)
#standard C+T
result.data.CT = result.data %>% 
  filter(r2_ind_vec==3&wc_ind_vec==1)
idx <- which.max(result.data.CT$r2.vec.test)
r2.max.ct <- result.data.CT$r2.vec.vad[idx]
idx <- which.max(r2.vec.test)
r2.max <- r2.vec.vad[idx]
r2.list.temp <- list(
  r2.max,
  r2.max.ct,
  result.data[idx,])
#save(r2.list.temp,file = paste0(out.dir,eth[i],"/r2.list_rho_eb_temp",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
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
tune = list(ntrees = c(10, 20, 50),
            max_depth = 1:3,
            shrinkage = c(0.001, 0.01, 0.1))

# Set detailed names = T so we can see the configuration for each function.
# Also shorten the name prefix.
learners = create.Learner("SL.xgboost", tune = tune, detailed_names = TRUE, name_prefix = "xgb")

SL.libray <- c(
  #"SL.xgboost"
  #"SL.randomForest"
  "SL.glmnet",
  "SL.ridge",
  #"SL.bayesglm"
  #"SL.stepAIC"
  "SL.nnet"
  #learners$names
  #,
  #"SL.svm"
  "SL.xgboost"
  #"SL.kernelKnn",
  #"SL.rpartPrune", 
  #"SL.lm"
  #"SL.mean"
)
options(mc.cores = 2)
getOption("mc.cores")

time1= proc.time()
sl = mcSuperLearner(Y = y.test, X = x.test, family = gaussian(),
                  # For a real analysis we would use V = 10.
                  # V = 3,
                  SL.library = SL.libray)
time2 = proc.time()-time1
print(time2)
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
result.data.CT = result.data %>% 
  filter(r2_ind_vec==3&wc_ind_vec==1)
idx <- which.max(result.data.CT$r2.vec.test)
r2.max.ct <- result.data.CT$r2.vec.vad[idx]
idx <- which.max(r2.vec.test)
r2.max <- r2.vec.vad[idx]
r2.list <- list(r2.stack,
                r2.max,
                r2.max.ct,
                result.data[idx,])
save(r2.list,file = paste0(out.dir,eth[i],"/r2.list_rho_eb_test_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
#write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump.auc.csv")
