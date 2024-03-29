#goal: rewrite simulation using plink2.0 to include both target and EUR
#use plink2 to calculate prs
#load LD_clump_file
#i is the ethnic
#l is the causal proportion
#m is the training sample size
#k is the p value thres
#j is the number of chromsome
args = commandArgs(trailingOnly = T)
i_rep = as.numeric(args[[1]])
i = as.numeric(args[[2]])
l = as.numeric(args[[3]])
m = as.numeric(args[[4]])
i1 = as.numeric(args[[5]])


eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
old.out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
#j = as.numeric(args[[3]])
sid <- Sys.getenv("SLURM_JOB_ID")
dir.create(paste0('/lscratch/',sid,'/test/'),showWarnings = F)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/prs/'),showWarnings = FALSE)
temp.dir.prs = paste0('/lscratch/',sid,'/test/prs/')

system(paste0("cp ",cur.dir,eth[i],"/all_chr_test.mega.* ",temp.dir))
system(paste0("ls ",temp.dir))
library(dplyr)
library(data.table)
setwd("/data/zhangh24/multi_ethnic/")

pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
#n.snp.mat <- matrix(0,length(pthres),4)


#for(l in 1:3){
print(l)
# for(m in 1:4){

print(m)
summary.eur <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[1],"/summary_out_rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)))  
colnames(summary.eur)[9] = "peur"
colnames(summary.eur)[7] = "beta_eur"
summary.eur.select = summary.eur %>% 
  mutate(sd_eur=beta_eur/STAT) %>% 
  select(SNP,A1,beta_eur,sd_eur,peur) %>% 
  rename(A1.EUR = A1)
r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)
sum.data <- as.data.frame(fread(paste0("./result/LD_simulation_GA/",eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))
colnames(sum.data)[2] <- "SNP"
#combine the target level summary stat with EUR
summary.com.match <- left_join(sum.data,summary.eur.select,by="SNP")

idx <- which(summary.com.match$A1!=summary.com.match$A1.EUR)
summary.com.match$A1.EUR[idx] <- summary.com.match$A1[idx]
summary.com.match$beta_eur[idx] <- -summary.com.match$beta_eur[idx]

summary.com.match = summary.com.match %>%
  mutate(beta_tar = BETA,
         sd_tar=BETA/STAT,
         z_stat_tar = STAT,
         z_stat_eur = beta_eur/sd_eur)
#estimate the prior
load(paste0(out.dir,eth[i],"/r2.list_rho_2DLD_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))

p.k1 =r2.list[[3]][[3]]
p.k2 = r2.list[[3]][[4]]
r_ind = r2.list[[3]][[5]]
w_ind = r2.list[[3]][[6]]
EstimatePrior <- function(beta_tar,sd_tar,
                          beta_eur,sd_eur){
  beta_tar = as.numeric(beta_tar)
  sd_tar = as.numeric(sd_tar)
  beta_eur = as.numeric(beta_eur)
  sd_eur = as.numeric(sd_eur)
  z_tar = beta_tar/sd_tar
  z_eur = beta_eur/sd_eur
  z_mat <-na.omit(cbind(z_tar,z_eur))
  
  
  prior.mat <- cov(z_mat)-diag(2)
  return(prior.mat)
}
EBpost <- function(beta_tar,sd_tar,
                   beta_eur,sd_eur,EBprior){
  
  beta_tar = as.numeric(beta_tar)
  sd_tar = as.numeric(sd_tar)
  beta_eur = as.numeric(beta_eur)
  sd_eur = as.numeric(sd_eur)
  prior.sigma = EBprior
  z_tar = beta_tar/sd_tar
  z_eur = beta_eur/sd_eur
  z_mat <-as.matrix(cbind(z_tar,z_eur))
  sd_mat =  as.matrix(cbind(sd_tar,sd_eur))
  post.sigma = solve(solve(prior.sigma)+diag(2))
  
  z_mat_post = z_mat
  
  p <- ncol(z_mat)
  
  for(k in 1:nrow(z_mat)){
    if(k%%10000==0){print(paste0(k," SNPs completed"))}
    z_temp = z_mat[k,]
    
    #find out nonmissing component
    
    idx <- which(!is.na(z_temp))
    if(length(idx)<p){
      z_temp <- z_temp[idx]
      
      post.sigma_temp = post.sigma[idx,idx,drop=F]
      z_post = post.sigma_temp%*%z_temp
    }else{
      z_post =post.sigma%*%z_temp
    }   
    
    z_mat_post[k,idx] = z_post
  }
  beta_mat_post = z_mat_post
  beta_mat_post[,1] =z_mat_post[,1]*sd_tar
  beta_mat_post[,2] =z_mat_post[,2]*sd_eur
  return(beta_mat_post)
}

LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
clump.snp <- LD
summary.com.prior = left_join(clump.snp,summary.com.match,by="SNP") %>% 
  filter(peur<p.k1|
           P<p.k2)
beta_tar <- summary.com.prior$beta_tar
sd_tar <- summary.com.prior$sd_tar
beta_eur <- summary.com.prior$beta_eur
sd_eur <- summary.com.prior$sd_eur

EBprior = EstimatePrior(beta_tar,sd_tar,
                        beta_eur,sd_eur)
beta_tar <- summary.com.match$beta_tar
sd_tar <- summary.com.match$sd_tar
beta_eur <- summary.com.match$beta_eur
sd_eur <- summary.com.match$sd_eur

post_beta_mat = EBpost(beta_tar,sd_tar,beta_eur,sd_eur,EBprior)
idx <- which(is.na(post_beta_mat[,2]))
post_beta_mat[idx,2] = 0
#post_beta_tar = post_beta_mat[,1,drop=F]



#summary.com.match$BETA = post_beta_tar


summary.com  = summary.com.match
summary.com = cbind(summary.com,post_beta_mat)
#remove duplicated snp
idx <- which(duplicated(summary.com$SNP))
if(length(idx)!=0){
  summary.com = summary.com[-idx,]
}
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    
    
    #read LD clumped SNPs
    LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    clump.snp <- LD[,1,drop=F] 
    #read the target ethnic group summary level statistics
    #combine the statistics with SNPs after clumping
    prs.all <- left_join(clump.snp,summary.com,by="SNP") 
    #idx <- which(summary.com$SNP=="rs201335322:154415264:C:T")
    colSums(is.na(prs.all))
    for(k1 in 1:length(pthres)){
      #keep al the SNPs with peur pass the threshold
      prs.all.temp = prs.all
      idx <- which(prs.all.temp$peur<=pthres[k1])
      prs.all.temp$P[idx] = 1E-20
      prs.file <- prs.all.temp %>% 
        select(SNP,A1,z_tar,z_eur,P)
      colSums(is.na(prs.file))
      write.table(prs.file,file = paste0(temp.dir.prs,"prs_file"),col.names = T,row.names = F,quote=F)
      
      #pdx <- which(duplicated(prs.file$SNP))
      
      
      
      p.value.file <- prs.all.temp %>% 
        select(SNP,P)
      write.table(p.value.file,file = paste0(temp.dir.prs,"p_value_file"),col.names = T,row.names = F,quote=F)
      n_pthres = length(pthres)
      q_range = data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres))
      temp = 1
      
      for(k2 in 1:length(pthres)){
        prs.file.sub <- prs.file %>% 
          filter(P<=pthres[k2]) %>% 
          select(SNP,A1,z_tar,z_eur,P)
        if(nrow(prs.file.sub)>0){
          q_range[temp,1] = paste0("p_value_",k2)
          q_range[temp,3] = pthres[k2]
          temp = temp+1
        }
      }
      q_range = q_range[1:(temp-1),]
      write.table(q_range,file = paste0(temp.dir.prs,"q_range_file"),row.names = F,col.names = F,quote=F)
      
      #remove NA
      res = system(paste0("/data/zhangh24/software/plink2_alpha --q-score-range ",
                          temp.dir.prs,"q_range_file ",temp.dir.prs,"p_value_file header --threads 2 ",
                          "--score-col-nums 3-4 ",
                          "--score ",temp.dir.prs,"prs_file cols=+scoresums,-scoreavgs ",
                          "--bfile ",temp.dir,"all_chr_test.mega ",
                          "--exclude ",old.out.dir,eth[i],"/duplicated.id  ",
                          "--out ",temp.dir.prs,"prs_eb_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1))
      print("step2 finished")
      #dup <- fread(paste0(old.out.dir,eth[i],"/duplicated.id"),header = F)
      
      
      #system(paste0("ls ",temp.dir.prs))
      if(res==3){
        stop()
      }
      #system(paste0("ls ",temp.dir.prs,"*.sscore"))
      res = system(paste0("mv ",temp.dir.prs,"*.sscore ",out.dir,eth[i],"/prs/"))
      if(res==2){
        stop()
      }
      system(paste0("rm -rf ",temp.dir.prs))
      dir.create(paste0(temp.dir.prs),showWarnings = FALSE)
      
      
      #system(paste0("/data/zhangh24/software/plink2 --score ",cur.dir,eth[i],"/prs/prs_file_pvalue_",k,"_rho_",l,"_size_",m,,"_rep_",i_rep," no-sum no-mean-imputation --bfile ",cur.dir,eth[i],"/all_chr.tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out ",cur.dir,eth[i],"/prs/prs_",k,"_rho_",l,"_size_",m))
    }
    
  }
  
}

# 
# n.test <- 10000
# n.vad <- n.test
# n.rep = 10
# #r2 mat represent the r2 matrix for the testing dataset
# #column represent the ethnic groups
# #row represent different p-value threshold
# cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
# setwd("/data/zhangh24/multi_ethnic/")
# out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
# y <- as.data.frame(fread(paste0(cur.dir,eth[i],"/phenotypes_rho",l,"_",i1,".phen")))
# y <- y[,2+(1:n.rep),drop=F]
# n <- nrow(y)
# y_test_mat <- y[(100000+1):nrow(y),,drop=F]
# 
# 
# r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
# wc_base_vec = c(50,100)
# 
# r2.vec.test <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
# r2.vec.vad <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
# pthres_vec1 <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
# pthres_vec2 <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
# r2_ind_vec <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
# wc_ind_vec <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
# prs.mat <- matrix(0,n.test+n.vad,2*length(pthres)^2*length(r2_vec)*length(wc_base_vec))
# temp = 1
# 
# for(r_ind in 1:length(r2_vec)){
#   wc_vec = wc_base_vec/r2_vec[r_ind]
#   for(w_ind in 1:length(wc_vec)){
#     print(c(r_ind,w_ind))
#     
#     LD <- as.data.frame(fread(paste0(out.dir,eth[i],"/LD_clump_two_way_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped")))
#     clump.snp <- LD
#     
#     #    colnames(sum.data)[2] <- "SNP"
#     
#     #for(k in 1:length(pthres)){
#     
#     prs.clump = left_join(clump.snp,summary.com,by="SNP")
#     
#     for(k1 in 1:length(pthres)){
#       for(k2 in 1:length(pthres)){
#         prs.all <- prs.clump %>%
#           filter(peur<=pthres[k1]|
#                    P<=pthres[k2])
#         if(nrow(prs.all)>0){
#           filename <- paste0(temp.dir.prs,"prs_eb_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,".p_value_",k2,".sscore")
#           
#           prs.temp <- fread(filename)  
#           prs.score <- prs.temp[,5:6]
#           prs.test <- prs.score[(1):(n.test)]
#           prs.vad <- prs.score[(n.test+1):(n.test+n.vad)]
#           prs.mat[,temp:(temp+1)] = as.matrix(prs.score)
#           temp = temp+2
#         }else{
#           prs.mat[,temp:(temp+1)] = 0
#           temp = temp+2
#         }
#         
#       }
#     }
#     
#   }
# }
# prs.sum = colSums(prs.mat)
# idx <- which(prs.sum!=0)
# #drop the prs with all 0
# prs.mat <- prs.mat[,idx]
# #drop the columns with perfect correlation
# prs.mat = as.data.frame(prs.mat)
# mtx = cor(prs.mat[1:n.test,])
# library(caret)
# drop = findCorrelation(mtx,cutoff=0.98)
# drop = names(prs.mat)[drop]
# prs.mat.new = prs.mat %>% 
#   select(-all_of(drop))
# library(SuperLearner)
# library(ranger)
# y.test = y_test_mat[1:n.test,i_rep]
# y.vad = y_test_mat[(n.test+1):(nrow(y_test_mat)),i_rep]
# x.test = as.data.frame(prs.mat.new[1:n.test,])
# x.vad= as.data.frame(prs.mat.new[(1+n.test):(n.test+n.vad),])
# 
# SL.libray <- c(
#   #"SL.xgboost"
#   #"SL.randomForest"
#   "SL.glmnet",
#   "SL.ridge",
#   #"SL.bayesglm"
#   #"SL.stepAIC"
#   "SL.nnet"
#   #"SL.ksvm",
#   #"SL.bartMachine", 
#   #"SL.kernelKnn",
#   #"SL.rpartPrune", 
#   #"SL.lm"
#   #"SL.mean"
# )
# sl = SuperLearner(Y = y.test, X = x.test, family = gaussian(),
#                   # For a real analysis we would use V = 10.
#                   # V = 3,
#                   SL.library = SL.libray)
# sl
# y.pred <- predict(sl, x.vad, onlySL = TRUE)
# #names(r2.vec.test) <- names(r2.vec.vad) <- pthres
# 
# #evaluate the best prs performance on the validation
# model <- lm(y.vad~y.pred[[1]])
# r2.stack <- summary(model)$r.square
# r2.list <- list(r2.stack)
# #update based on including both eur and target coefficients
# save(r2.list,file = paste0(out.dir,eth[i],"/r2.list_update_rho_eb_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
# 
# #pthres <- c(1E-10,1E-09,5E-08,1E-07,2.5E-07,5E-07,7.5E-07,1E-06,2.5E-06,5E-06,7.5E-06,1E-05,2.5e-05,5E-05,7.5e-05,1E-04,2.5E-04,5E-04,7.5E-04,1E-03)
# 
# system(paste0("rm -rf ", temp.dir))
# 
# #}
# 
# 
# 
# a = rep(0,10)
# for(i_rep in 1:10){
#   load(paste0(out.dir,eth[i],"/r2.list_rho_eb_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))
#   a[i_rep] = r2.list[[1]]
# }




