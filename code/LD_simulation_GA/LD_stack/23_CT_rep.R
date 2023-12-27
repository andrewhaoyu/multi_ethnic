#calculate AUC for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i1 = as.numeric(args[[4]])
library(boot)
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

r2_vec = c(0.1)
wc_base_vec = c(50)
n.rep = 10
ci_low = rep(0,10)
ci_high = rep(0,10)
r2_vad = rep(0,10)
r2.vec.test <- rep(0,length(pthres)*length(r2_vec)*length(wc_base_vec)*n.rep)
r2.vec.vad <- rep(0,length(pthres)*length(r2_vec)*length(wc_base_vec)*n.rep)
pthres_vec <- rep(0,length(pthres)*length(r2_vec)*length(wc_base_vec)*n.rep)
r2_ind_vec <- rep(0,length(pthres)*length(r2_vec)*length(wc_base_vec)*n.rep)
wc_ind_vec <- rep(0,length(pthres)*length(r2_vec)*length(wc_base_vec)*n.rep)
rep_vec = rep(0,length(pthres)*length(r2_vec)*length(wc_base_vec)*n.rep)
temp = 1
for(i_rep in 1:n.rep){
  sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
  
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
          
          
          filename <- paste0(out.dir,eth[i],"/prs/prs_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".p_value_",k,".profile")
          
          prs.temp <- fread(filename)
          prs.score <- prs.temp$SCORE
          #prs.score <- prs.temp$SCORE*2*length(idx)+prs.score
          #prs.score.mat[,k] = prs.score
          prs.test <- prs.score[(1):(n.test)]
          prs.vad <- prs.score[(n.test+1):(n.test+n.vad)]
          #model = lm(y~prs.score)
          y.test = y_test_mat[1:n.test,i_rep]
          
          model1 <- lm(y.test~prs.test)
          #r2.test.rep[i_rep] <- summary(model1)$r.square
          
          r2.vec.test[temp] = summary(model1)$r.square
          
          pthres_vec[temp] = pthres[k]
          r2_ind_vec[temp] = r_ind
          wc_ind_vec[temp] = w_ind
          rep_vec[temp] = i_rep
          
          temp = temp+1
        }else{
          r2.vec.test[temp] = 0
          pthres_vec[temp] = pthres[k]
          r2_ind_vec[temp] = r_ind
          wc_ind_vec[temp] = w_ind
          rep_vec[temp] = i_rep
          temp = temp+1
        }
        
      }
    }
    
  }
  r2_function <- function(data, indices) {
    sample_data <- data[indices, ]
    model <- lm(y.vad ~ prs.vad, data = sample_data)
    return(summary(model)$r.square)
  }
  #find best CT cutoff
  idx = which.max(r2.vec.test)
  prs.all <- prs.clump %>% 
    filter(P<=pthres[idx])
  filename <- paste0(out.dir,eth[i],"/prs/prs_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".p_value_",k,".profile")
  
  prs.temp <- fread(filename)
  prs.vad <- prs.score[(n.test+1):(n.test+n.vad)]
  y.vad = y_test_mat[(n.test+1):(nrow(y_test_mat)),i_rep]
  model2 <- lm(y.vad~prs.vad)
  data <- data.frame(y.vad = y.vad, prs.vad = prs.vad)  # assuming y.vad and prs.vad are your data vectors
  # results <- boot(data = data, statistic = r2_function, R = 1000)
  # ci_result = boot.ci(results, type = "perc")
  # ci_low[i_rep] = ci_result$percent[4]
  # ci_high[i_rep] = ci_result$percent[5]
  r2_vad[i_rep] = summary(model2)$r.square
}




#save the analysis table with different p-value threshold
result.data = data.frame(r2_vad)



save(result.data,file = paste0(out.dir,eth[i],"/ctrep_",l,"_size_",m,"_GA_",i1))
#}
#}
#write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump.auc.csv")
