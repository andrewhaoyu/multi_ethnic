args = commandArgs(trailingOnly = T)
t_rep = as.numeric(args[[1]])
#test the running time and memory for TDLD-SLEB and PRS-CSx on CHR 22
time_ct_sleb_start = proc.time()
i = 2
library(data.table)
library(dplyr)
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/LD/'),showWarnings = FALSE)
temp.dir.LD <- paste0('/lscratch/',sid,'/test/LD/')
eth <- c("EUR","AFR","AMR","EAS","SAS")
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/time_mem/"
library(CTSLEB)
#load EUR summary data
load(paste0(out.dir,"eur_sumdata.rdata"))
#load AFR summary data
load(paste0(out.dir,"tar_sumdata.rdata"))
summary.eur = summary.eur %>% 
  mutate(SE = BETA/STAT) %>% 
  rename(P = peur)
summary.tar = summary.tar %>% 
  mutate(SE = BETA/STAT)

sum_com <- AlignSum(sum_tar = summary.tar,
                    sum_other = summary.eur)
#split the SNPs into two groups
sum_com_split <- SplitSum(sum_com)
#sum_com_split is a list with two data frame
#the first data frame contains SNPs with p_eur < p_target, the P value column is from p_eur
sum_other_ref = sum_com_split[[1]]
#the second data frame contains target population-specific SNPs or p_eur < p_target. The p-value column is from p_target. 
sum_tar_ref = sum_com_split[[2]]

r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)
data.dir = temp.dir 
soft.dir = "/data/zhangh24/software/"
write.table(sum_other_ref,paste0(temp.dir,"sum_other_ref"),col.names = T,row.names = F,quote=F)
write.table(sum_tar_ref,paste0(temp.dir,"sum_tar_ref"),col.names = T,row.names = F,quote=F)
system(paste0("cp ",out.dir,"EUR_ref_chr22.bed ",temp.dir,"EUR_ref_chr22.bed"))
system(paste0("cp ",out.dir,"EUR_ref_chr22.bim ",temp.dir,"EUR_ref_chr22.bim"))
system(paste0("cp ",out.dir,"EUR_ref_chr22.fam ",temp.dir,"EUR_ref_chr22.fam"))
system(paste0("cp ",out.dir,"AFR_ref_chr22.bed ",temp.dir,"AFR_ref_chr22.bed"))
system(paste0("cp ",out.dir,"AFR_ref_chr22.bim ",temp.dir,"AFR_ref_chr22.bim"))
system(paste0("cp ",out.dir,"AFR_ref_chr22.fam ",temp.dir,"AFR_ref_chr22.fam"))

setwd("/data/zhangh24/multi_ethnic/")
snp_list = list()
temp = 1
for(r_ind in 1:length(r2_vec)){
  #create the window size given the clumping r2
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    pthr = 1
    r2thr = r2_vec[r_ind]
    kbpthr = wc_vec[w_ind]
    #for the first group, we perform clumping using EUR as reference   
    system(paste0(soft.dir,"plink2 ",
                  "--bfile ",data.dir,"EUR_ref_chr22 ",
                  "--clump ",temp.dir,"sum_other_ref ",
                  "--clump-p1 ",pthr," ",
                  "--clump-r2 ",r2thr," ",
                  "--clump-kb ",kbpthr," ",
                  "--out ", temp.dir,"EUR_ref_CT_rind_",r_ind,"_wcind_",w_ind))
    #for the second group, we perform clumping using AFR as reference
    system(paste0(soft.dir,"plink2 ",
                  "--bfile ",data.dir,"AFR_ref_chr22 ",
                  "--clump ",temp.dir,"sum_tar_ref ",
                  "--clump-p1 ",pthr," ",
                  "--clump-r2 ",r2thr," ",
                  "--clump-kb ",kbpthr," ",
                  "--out ", temp.dir,"AFR_ref_CT_rind_",r_ind,"_wcind_",w_ind))
    
    #combine the SNPs from the two clumping groups
    LD_EUR= fread(paste0(temp.dir,"EUR_ref_CT_rind_",r_ind,"_wcind_",w_ind,".clumped"))[,3,drop=F]
    LD_tar = fread(paste0(temp.dir,"AFR_ref_CT_rind_",r_ind,"_wcind_",w_ind,".clumped"))[,3,drop=F]
    LD  = rbind(LD_EUR,LD_tar)
    snp_list[[temp]] = LD
    names(snp_list[[temp]]) = paste0("clump_r2_",r2thr,"_ws_",kbpthr)
    temp = temp + 1
  }
}
time_ct_end = proc.time()
#combine the LD clumping results
#calculate prs
#create q-value file for prs analysis
system(paste0("cp ",out.dir,"AFR_test_mega_chr22.bed ",temp.dir,"AFR_test_mega_chr22.bed"))
system(paste0("cp ",out.dir,"AFR_test_mega_chr22.bim ",temp.dir,"AFR_test_mega_chr22.bim"))
system(paste0("cp ",out.dir,"AFR_test_mega_chr22.fam ",temp.dir,"AFR_test_mega_chr22.fam"))


plink_file = PreparePlinkFile(snp_list,sum_com)

#score_file description
#the first column contains the unique SNPs after clumping results under all combinations of r2-cutoff and window_size
#the second column is the effect allele
#the third to the last columns contains the regression coefficients of the target population for SNPs after LD-clumping under a specific combination of r2-cutoff and base_window_size
#the coefficients is put as 0 if a SNP doesn't exist in the clumping results under a specific combination of r2-cutoff and base_window_size
score_file = plink_file[[1]]
write.table(score_file,file = paste0(temp.dir,"score_file"),row.names = F,col.names = F,quote=F)
#p_value_file description
#the first column is the same as score_file
#the second column is the p-values of SNPs from the GWAS of the target population
p_value_file = plink_file[[2]]
# unique_infor description
#unique_infor contains the information for all SNPs after the clumping step
#unique_infor has SNP, CHR, BP, A1 (effect_allele), 
#(BETA, SE, P) for the target population, 
#(BETA_other, SE_other,P_other) for the EUR  population
unique_infor = plink_file[[3]]
#specific p-value threshold
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
#create q-range file
q_range = CreateQRange(pthres)

write.table(q_range,file = paste0(temp.dir,"q_range_file"),row.names = F,col.names = F,quote=F)

#create a temporary p_value_file
p_value_file_temp = p_value_file
for(k1 in 1:length(pthres)){
  #keep al the SNPs with P_EUR less than pthres[k1] in the analyses
  idx <- which(unique_infor$P_other<=pthres[k1])
  p_value_file_temp$P[idx] = 0
  write.table(p_value_file_temp,file = paste0(temp.dir,"p_value_file"),col.names = F,row.names = F,quote=F)
  n_col = ncol(score_file)
  #the output of plink2 create 9 different files named as prs_p_other_k1.p_tar_k2.sscore
  #this output file contains 16 columns
  #the column 1-4 are: family ID, individual ID, 2*total number of SNPs in the PRS, the sum of allele count
  #column 5-16 are the PRS scores with SNP of p_target<p_thres[k2]|p_eur<p_thres[k1] for different combinations of r2-cutoff and base_window_size
  res = system(paste0(soft.dir,"plink2_alpha ",
                      "--q-score-range ",temp.dir,"q_range_file ",temp.dir,"p_value_file ",
                      "--score-col-nums 3-",n_col," ",
                      "--score ",temp.dir,"score_file  ",
                      "--bfile ",data.dir,"AFR_test_chr22 ",
                      "--out ",temp.dir,"prs_p_other_",k1))
  
}
#combine all the prs
prs_list = list()
temp = 1
#take the column name of different clumping parameters
names = colnames(score_file)[3:ncol(score_file)]
for(k1 in 1:length(pthres)){
  for(k2 in 1:length(pthres)){
    #the --score command in plink2 computes PRS as G*beta/(2*number of SNPs)
    #This doesn't affect prediction R2 of the PRS if you compute PRS for all chromosomes together
    #But if you compute PRS by chromosome, you need to times (2*number of SNPs) before sum the score together
    #here I am scaling the score by (2*number of SNPs) to aviod this potential
    
    #load PRS for SNPs with p_target<p_thres[k2]|p_eur<p_thres[k1] 
    prs_temp = fread(paste0(temp.dir,"prs_p_other_",k1,".p_tar_",k2,".sscore"))
    # times (2*number of SNPs)
    prs_list[[temp]] = prs_temp[,5:ncol(prs_temp)]*2*nrow(score_file)
    
    colnames(prs_list[[temp]]) = paste0(names,"_","p_other_",pthres[k1],"_p_tar_",pthres[k2])
    temp = temp + 1
  }
}
prs_mat = as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
prs_tun = prs_mat[1:10000,]

#find the best R-square among the all the PRSs to find candidate set
#we use this candidate for estimating covariance matrix for the prior distribution
#create prediction r2 vector to store r2 for different prs
n.total.prs = length(pthres)^2*length(r2_vec)*length(wc_base_vec)
prs_r2_vec_test = rep(0,n.total.prs)
time_prs1_end = proc.time()


#find the optimal r2 performance based on testing data
load(paste0(out.dir,"y_test.rdata"))
for(p_ind in 1:n.total.prs){
  #the first two columns of prs_tun are family id and individual id
  #prs starts from the third column
  model = lm(y_test~prs_tun[,(2+p_ind)])
  prs_r2_vec_test[p_ind] = summary(model)$r.square
}
max_ind = which.max(prs_r2_vec_test)
#+2 is due to the first two columns are family id and individual id
print(colnames(prs_tun)[max_ind+2])
r2.vec.test <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
pthres_vec1 <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
pthres_vec2 <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
r2_ind_vec <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
wc_ind_vec <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))

temp = 1
n.test = 10000
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    for(k1 in 1:length(pthres)){
      for(k2 in 1:length(pthres)){
        
        filename <- paste0(temp.dir.prs,"prs_2DLD_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,".p_value_",k2,".profile")
        prs.temp <- fread(filename)  
        prs.score <- prs.temp$SCORE
        prs.test <- prs.score[(1):(n.test)]
        model1 <- lm(y_test~prs.test)
        r2.vec.test[temp] = summary(model1)$r.square
        pthres_vec1[temp] = pthres[k1]
        pthres_vec2[temp] = pthres[k2]
        r2_ind_vec[temp] = r_ind
        wc_ind_vec[temp] = w_ind
        temp = temp + 1
      }
    }
  }
}
idx <- which.max(r2.vec.test)
p.k1 =pthres_vec1[idx]
p.k2 = pthres_vec2[idx]
r_ind = r2_ind_vec[idx]
w_ind = wc_ind_vec[idx]
time_first_compute_r2_end = proc.time()
#compute the prs based on TDLD clumping results
source("/data/zhangh24/multi_ethnic/code/stratch/EB_function.R")
LD.EUR <- as.data.frame(fread(paste0(temp.dir,eth[1],"_LD_clump_two_dim_rind_",r_ind,"_wcind_",w_ind,".clumped")))
LD.tar <- as.data.frame(fread(paste0(temp.dir,eth[i],"_LD_clump_two_dim_rind_",r_ind,"_wcind_",w_ind,".clumped")))

LD.EUR <- LD.EUR[,3,drop=F]
LD.tar <- LD.tar[,3,drop=F]
LD <- rbind(LD.EUR,LD.tar)

#combine the statistics with SNPs after clumping
load(paste0(out.dir,"other_ethnic_beta_mat"))
load(paste0(out.dir,"other_ethnic_sd_mat"))
summary.com$beta_mat = beta.mat
summary.com$sd_mat = sd.mat
summary.com.prior = left_join(LD,summary.com,by="SNP") %>% 
  filter(peur<p.k1|
           P<p.k2)

beta_tar = summary.com.prior$beta_tar
sd_tar = summary.com.prior$sd_tar
beta_eur = summary.com.prior$beta_eur
sd_eur = summary.com.prior$sd_eur
beta_other_mat = summary.com.prior$beta_mat
sd_other_mat = summary.com.prior$sd_mat

prior.sigma = EstimatePriorMulti(beta_tar,sd_tar,
                                 beta_eur,sd_eur,
                                 beta_other_mat,
                                 sd_other_mat)

beta_tar <- summary.com$beta_tar
sd_tar <- summary.com$sd_tar
beta_eur <- summary.com$beta_eur
sd_eur <- summary.com$sd_eur
beta_other_mat <- summary.com$beta_mat
sd_other_mat <- summary.com$sd_mat

post_beta_tar = EBpostMulti(beta_tar,sd_tar,
                            beta_eur,sd_eur,beta_other_mat,
                            sd_other_mat,
                            prior.sigma)[,1,drop=F]
summary.com$BETA = post_beta_tar
time_eb_end = proc.time()
#calculate the prs using TDLD-EB
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    LD.EUR <- as.data.frame(fread(paste0(temp.dir,eth[1],"_LD_clump_two_dim_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    LD.tar <- as.data.frame(fread(paste0(temp.dir,eth[i],"_LD_clump_two_dim_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    
    LD.EUR <- LD.EUR[,3,drop=F]
    LD.tar <- LD.tar[,3,drop=F]
    LD <- rbind(LD.EUR,LD.tar)
    
    #combine the statistics with SNPs after clumping
    prs.all <- left_join(LD,summary.com,by="SNP") 
    colSums(is.na(prs.all))
    
    for(k1 in 1:length(pthres)){
      #keep al the SNPs with peur pass the threshold
      prs.file = prs.all %>% 
        mutate(P = replace(P,peur<=pthres[k1],1E-20)) %>% 
        select(SNP,A1,BETA,P)
      write.table(prs.file,file = paste0(temp.dir.prs,"prs_file"),col.names = T,row.names = F,quote=F)
      
      p.value.file <- prs.file %>% 
        select(SNP,P)
      write.table(p.value.file,file = paste0(temp.dir.prs,"p_value_file"),col.names = T,row.names = F,quote=F)
      
      old.out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
      res = system(paste0("/data/zhangh24/software/plink2 --q-score-range ",temp.dir.prs,"q_range_file ",temp.dir.prs,"p_value_file header --threads 2 --score ",temp.dir.prs,"prs_file header no-sum no-mean-imputation --bfile ",temp.dir,"AFR_test_mega_chr22 --exclude ",old.out.dir,eth[i],"/duplicated.id --out ",temp.dir.prs,"prs_eb_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1))
      print("step2 finished")
      
      
      #system(paste0("/data/zhangh24/software/plink2 --score ",cur.dir,eth[i],"/prs/prs_file_pvalue_",k,"_rho_",l,"_size_",m,,"_rep_",i_rep," no-sum no-mean-imputation --bfile ",cur.dir,eth[i],"/all_chr.tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out ",cur.dir,eth[i],"/prs/prs_",k,"_rho_",l,"_size_",m))
    }
    
  }
  
  
}

load(paste0(out.dir,"y_test.rdata"))
r2.vec.test <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
pthres_vec1 <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
pthres_vec2 <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
r2_ind_vec <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
wc_ind_vec <- rep(0,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
prs.mat <- matrix(0,n.test,length(pthres)^2*length(r2_vec)*length(wc_base_vec))
temp = 1
n.test = 10000
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    for(k1 in 1:length(pthres)){
      for(k2 in 1:length(pthres)){
        
        filename <- paste0(temp.dir.prs,"prs_2DLD_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1,".p_value_",k2,".profile")
        prs.temp <- fread(filename)  
        prs.score <- prs.temp$SCORE
        prs.test <- prs.score[(1):(n.test)]
        model1 <- lm(y_test~prs.test)
        r2.vec.test[temp] = summary(model1)$r.square
        pthres_vec1[temp] = pthres[k1]
        pthres_vec2[temp] = pthres[k2]
        r2_ind_vec[temp] = r_ind
        wc_ind_vec[temp] = w_ind
        prs.mat[,temp] = prs.test
        temp = temp + 1
      }
    }
  }
}
time_prs2_end = proc.time()
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
prs.mat = prs.mat %>% 
  select(-all_of(drop))

library(SuperLearner)
library(ranger)
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
sl = SuperLearner(Y = y_test, X = prs.mat, family = gaussian(),
                  # For a real analysis we would use V = 10.
                  # V = 3,
                  SL.library = SL.libray)
#y.pred <- predict(sl, x.vad, onlySL = TRUE)
time_super_learning_end = proc.time()
y.pred <- predict(sl, prs.mat, onlySL = TRUE)
model = lm(y_test~y.pred[[1]])
time_ct_sleb_end = proc.time()
total_time = time_ct_sleb_end-time_ct_sleb_start
train_time = time_ct_end - time_ct_sleb_start + time_eb_end - time_first_compute_r2_end + time_super_learning_end - time_prs2_end
prs_time = time_prs1_end - time_ct_end + time_prs2_end - time_eb_end
compute_r2_time = time_first_compute_r2_end - time_prs1_end + time_ct_sleb_end - time_super_learning_end
time_vec = rbind(total_time,train_time,prs_time,compute_r2_time)
save(time_vec,file = paste0(out.dir,"TDLD_SLEBalleth_trep_",t_rep,".rdata"))
