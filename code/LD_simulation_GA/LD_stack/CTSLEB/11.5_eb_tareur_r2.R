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

prs.mat <- matrix(0,n.test+n.vad,2*length(pthres)^2*length(r2_vec)*length(wc_base_vec))
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
          filename <- paste0(out.dir,eth[i],"/prs/prs_eb_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,".p_value_",k2,".sscore")

          prs.temp <- fread(filename)
          prs.score <- prs.temp[,5:6]
          prs.test <- prs.score[(1):(n.test)]
          prs.vad <- prs.score[(n.test+1):(n.test+n.vad)]
          prs.mat[,temp:(temp+1)] = as.matrix(prs.score)
          temp = temp+2
        }else{
          prs.mat[,temp:(temp+1)] = 0
          temp = temp+2
        }

      }
    }

  }
}
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
y.test = y_test_mat[1:n.test,i_rep]
y.vad = y_test_mat[(n.test+1):(nrow(y_test_mat)),i_rep]
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
sl
y.pred <- predict(sl, x.vad, onlySL = TRUE)
#names(r2.vec.test) <- names(r2.vec.vad) <- pthres

#evaluate the best prs performance on the validation
model <- lm(y.vad~y.pred[[1]])
r2.stack <- summary(model)$r.square
r2.list <- list(r2.stack)
#update based on including both eur and target coefficients
save(r2.list,file = paste0(out.dir,eth[i],"/r2.list_update_rho_eb_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
