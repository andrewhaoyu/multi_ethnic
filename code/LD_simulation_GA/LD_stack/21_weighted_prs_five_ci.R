#calculate the AUC for EUR SNPs with EUR coefficients
#with target ethnic coefficients
#with EB coefficients.
#calculate AUC for different ethnic groups
#i for different ethnic groups
#l for different causal proportion
#m for different traning sample size
#q for three different methods

#i_rep = as.numeric(args[[5]])
args = commandArgs(trailingOnly = T)

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
i_rep = as.numeric(args[[3]])
#j = as.numeric(args[[4]])
#l = as.numeric(args[[3]])
m = as.numeric(args[[4]])
i1 = as.numeric(args[[5]])
eth <- c("EUR","AFR","AMR","EAS","SAS")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)

cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
sid <- Sys.getenv("SLURM_JOB_ID")
dir.create(paste0('/lscratch/',sid,'/test/'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/prs/'),showWarnings = FALSE)
temp.dir.prs = paste0('/lscratch/',sid,'/test/prs/')
system(paste0("cp ",cur.dir,eth[i],"/all_chr_test.mega.* ",temp.dir,"."))
system(paste0("ls ",temp.dir))
method <- c("Weighted-PRS (five ancestries)")

library(dplyr)
library(data.table)
#n <- 120000
n.train.vec <- c(15000,45000,80000,100000)
n.rep = 10

#r2 mat represent the r2 matrix for the testing dataset
#column represent the ethnic groups
#row represent different p-value threshold

pheno.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
setwd("/data/zhangh24/multi_ethnic/")
n.test = 10000
r_ind = 1
w_ind = 1
#create a matrix to save PRS
prs_com = matrix(0,20000,5)
#prepare the coefficients for best CT PRS in each non-EUR ans
for(i_eth in 2:5){
  y <- as.data.frame(fread(paste0(pheno.dir,eth[i_eth],"/phenotypes_rho",l,"_",i1,".phen")))
  y <- y[,2+(1:n.rep),drop=F]
  n <- nrow(y)
  y_test_mat <- y[(100000+1):nrow(y),,drop=F]
  y.test = y_test_mat[1:n.test,i_rep]
  sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i_eth],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))  
  LD <- as.data.frame(fread(paste0(out.dir,eth[i_eth],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
  clump.snp <- LD[,3,drop=F]  
  
  
  
  colnames(sum.data)[2] <- "SNP"
  #for(k in 1:length(pthres)){
  
  #for(k in 1:length(pthres)){
  
  prs.clump = left_join(clump.snp,sum.data,by="SNP")
  temp = 1
  r2.vec.test = rep(0, length(pthres))
  for(k in 1:length(pthres)){
    prs.all <- prs.clump %>% 
      filter(P<=pthres[k])
    if(nrow(prs.all)>0){
      
      
      filename <- paste0(out.dir,eth[i_eth],"/prs/prs_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".p_value_",k,".profile")
      
      prs.temp <- fread(filename)
      prs.score <- prs.temp$SCORE
      #prs.score <- prs.temp$SCORE*2*length(idx)+prs.score
      #prs.score.mat[,k] = prs.score
      prs.test <- prs.score[(1):(n.test)]
      model1 <- lm(y.test~prs.test)
      #r2.test.rep[i_rep] <- summary(model1)$r.square
      
      r2.vec.test[temp] = summary(model1)$r.square
      temp = temp+1
    }else{
      r2.vec.test[temp] = 0
    
      temp = temp+1
    }
  }
  #find the best cutoff
  idx = which.max(r2.vec.test)
  #calculate the prs
  prs.coef = prs.clump %>% 
    filter(P<=pthres[idx]) %>% select(SNP,A1,BETA)
  write.table(prs.coef,file = paste0(temp.dir.prs,"prs_file"),col.names = T,row.names = F,quote=F)
  res <- system(paste0("/data/zhangh24/software/plink2 ",
                       "--threads 2 --score ",temp.dir.prs,"prs_file header no-sum no-mean-imputation ",
                       "--bfile ",temp.dir,"all_chr_test.mega ",
                       "--exclude /data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/duplicated.id  ",
                       "--out ",temp.dir.prs,"prs_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind))
  prs_score = fread(paste0(temp.dir.prs,"prs_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".profile"))
  prs_com[,i_eth] = prs_score$SCORE
}
#load EUR PRS
filename <- paste0(out.dir,eth[i],"/prs/prs_besteur_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".sscore")
prs.temp <- as.data.frame(fread(filename))
#load best EUR PRS
q = 1
prs.score.eur<- prs.temp[,4+q]
prs_com[,1] = prs.score.eur

y <- as.data.frame(fread(paste0(pheno.dir,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
y <- y[,2+(1:n.rep),drop=F]
n <- nrow(y)
y_test_mat <- y[(100000+1):nrow(y),,drop=F]
n.vad = 10000
y.test = y_test_mat[1:n.test,i_rep]
y.vad = y_test_mat[(1+n.test):(n.test+n.vad),i_rep]
prs.score.test = prs_com[1:n.test,]
prs.score.vad = prs_com[(1+n.test):(n.test+n.vad),]
model1 <- lm(y.test~prs.score.test)

beta = coefficients(model1)[2:6]
model2 = lm(y.vad~prs.score.vad%*%beta)
library(boot)
r2 = summary(model2)$r.square
r2_function <- function(data, indices) {
  sample_data <- data[indices, ]
  model <- lm(y.vad ~ prs.vad, data = sample_data)
  return(summary(model)$r.square)
}
data <- data.frame(y.vad = y.vad, prs.vad = prs.score.vad%*%beta)  # assuming y.vad and prs.vad are your data vectors
results <- boot(data = data, statistic = r2_function, R = 1000)
ci_result = boot.ci(results, type = "perc")
ci_low= ci_result$percent[4]
ci_high = ci_result$percent[5]

weightedprs.result = data.frame(r2 = r2,
                                ci_low = ci_low,
                                ci_high = ci_high)

save(weightedprs.result,file = paste0(out.dir,
                                      "weighted_prs_five_ci_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".rdata"))







#load(paste0(out.dir,"LD.clump.result.GA_",i1,".rdata"))


#save(r2.vec,file = paste0(out.dir,eth[i],"/r2_eursnp_",method[q],"_rho_",l,"_size_",m,"_GA_",i1))
#write.csv(r2.mat,file = "/data/zhangh24/multi_ethnic/result/LD_simulation/ld.clump.auc.csv")
