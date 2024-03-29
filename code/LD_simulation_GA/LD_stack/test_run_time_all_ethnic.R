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
#load EUR summary data
load(paste0(out.dir,"eur_sumdata.rdata"))
#only keep the id and p-value for mathching
summary.eur.select = summary.eur %>% 
  rename(beta_eur = BETA) %>% 
  mutate(sd_eur=beta_eur/STAT) %>% 
  select(SNP,A1,beta_eur,sd_eur,peur) %>% 
  rename(A1.EUR = A1)
#load AFR summary data
load(paste0(out.dir,"tar_sumdata.rdata"))
summary.com <- left_join(summary.tar,summary.eur.select,by="SNP")
#match alleles
idx <- which(summary.com$A1!=summary.com$A1.EUR)
summary.com$A1.EUR[idx] <- summary.com$A1[idx]
summary.com$beta_eur[idx] <- -summary.com$beta_eur[idx]
summary.com = summary.com %>% 
  rename(beta_tar = BETA) %>% 
  mutate(sd_tar = beta_tar/STAT)
#select the SNPs from EUR p-value
idx <- which(summary.com$peur<summary.com$P)
#generate assoc files for clumping using EUR as reference
assoc =  summary.com[idx,] %>%
  select(SNP,peur) %>% 
  rename(P=peur
  )

write.table(assoc,paste0("/lscratch/",sid,"/test/",eth[1],"_assoc.out"),col.names = T,row.names = F,quote=F)

system(paste0("cp ",out.dir,"EUR_ref_chr22.bed ",temp.dir,"EUR_ref_chr22.bed"))
system(paste0("cp ",out.dir,"EUR_ref_chr22.bim ",temp.dir,"EUR_ref_chr22.bim"))
system(paste0("cp ",out.dir,"EUR_ref_chr22.fam ",temp.dir,"EUR_ref_chr22.fam"))
setwd("/data/zhangh24/multi_ethnic/")
r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)
n.test = 10000
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    pthr = 0.5
    r2thr = r2_vec[r_ind]
    kbpthr = wc_vec[w_ind]
    eth <- c("EUR","AFR","AMR","EAS","SAS")
    
    #code <- rep("c",5*3*3)
    #system(paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/",eth[i],"/chr_all --clump ",cur.dir,eth[i],"/summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,".out --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep))
    res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",temp.dir,"EUR_ref_chr22 --clump ",temp.dir,eth[1],"_assoc.out --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ", temp.dir,eth[1],"_LD_clump_two_dim_rind_",r_ind,"_wcind_",w_ind))
    if(res==2){
      stop()
    }
    #system(paste0("mv ",temp.dir,"LD_clump_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,"_rind_",r_ind,"_wcind_",w_ind,".clumped ",out.dir,eth[i],"/"))
  }
}

system(paste0("rm ",temp.dir,"*.bim"))
system(paste0("rm ",temp.dir,"*.bed"))
system(paste0("rm ",temp.dir,"*.fam"))

idx <- which((summary.com$P<summary.com$peur)|is.na(summary.com$peur))
assoc <- summary.com[idx,] %>% 
  select(SNP,P) 
#%>% 
#rename(P=peur)
write.table(assoc,paste0("/lscratch/",sid,"/test/",eth[i],"_assoc.out"),col.names = T,row.names = F,quote=F)
system(paste0("cp ",out.dir,"AFR_ref_chr22.bed ",temp.dir,"AFR_ref_chr22.bed"))
system(paste0("cp ",out.dir,"AFR_ref_chr22.bim ",temp.dir,"AFR_ref_chr22.bim"))
system(paste0("cp ",out.dir,"AFR_ref_chr22.fam ",temp.dir,"AFR_ref_chr22.fam"))


for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    print(c(r_ind,w_ind))
    pthr = 0.5
    r2thr = r2_vec[r_ind]
    kbpthr = wc_vec[w_ind]
    eth <- c("EUR","AFR","AMR","EAS","SAS")
    
    res = system(paste0("/data/zhangh24/software/plink2 --threads 2 ",
                        "--bfile ",temp.dir,"AFR_ref_chr22 ",
                        "--clump ",temp.dir,eth[i],"_assoc.out ",
                        "--clump-p1 ",pthr," --clump-r2 ",r2thr,"  ",
                        "--clump-kb ",kbpthr," --out ", temp.dir,eth[i],"_LD_clump_two_dim_rind_",r_ind,"_wcind_",w_ind))
    if(res==2){
      stop()
    }
    
  }
}
system(paste0("rm ",temp.dir,"*.bim"))
system(paste0("rm ",temp.dir,"*.bed"))
system(paste0("rm ",temp.dir,"*.fam"))
time_ct_end = proc.time()
dir.create(paste0('/lscratch/',sid,'/test/prs/'),showWarnings = FALSE)
temp.dir.prs = paste0('/lscratch/',sid,'/test/prs/')
#combine the LD clumping results
#calculate prs
#create q-value file for prs analysis
system(paste0("cp ",out.dir,"AFR_test_mega_chr22.bed ",temp.dir,"AFR_test_mega_chr22.bed"))
system(paste0("cp ",out.dir,"AFR_test_mega_chr22.bim ",temp.dir,"AFR_test_mega_chr22.bim"))
system(paste0("cp ",out.dir,"AFR_test_mega_chr22.fam ",temp.dir,"AFR_test_mega_chr22.fam"))


pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01)
n_pthres = length(pthres)
q_range = data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres))
temp = 1

for(k2 in 1:length(pthres)){
  
  
  q_range[temp,1] = paste0("p_value_",k2)
  q_range[temp,3] = pthres[k2]
  temp = temp+1
  
}
q_range = q_range[1:(temp-1),]
write.table(q_range,file = paste0(temp.dir.prs,"q_range_file"),row.names = F,col.names = F,quote=F)


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
        rename(BETA=beta_tar) %>% 
        select(SNP,A1,BETA,P)
      write.table(prs.file,file = paste0(temp.dir.prs,"prs_file"),col.names = T,row.names = F,quote=F)
      
      p.value.file <- prs.file %>% 
        select(SNP,P)
      write.table(p.value.file,file = paste0(temp.dir.prs,"p_value_file"),col.names = T,row.names = F,quote=F)
      
      old.out.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
      res = system(paste0("/data/zhangh24/software/plink2 --q-score-range ",temp.dir.prs,"q_range_file ",temp.dir.prs,"p_value_file header --threads 2 --score ",temp.dir.prs,"prs_file header no-sum no-mean-imputation --bfile ",temp.dir,"AFR_test_mega_chr22 --exclude ",old.out.dir,eth[i],"/duplicated.id --out ",temp.dir.prs,"prs_2DLD_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1))
      print("step2 finished")
      
      
      #system(paste0("/data/zhangh24/software/plink2 --score ",cur.dir,eth[i],"/prs/prs_file_pvalue_",k,"_rho_",l,"_size_",m,,"_rep_",i_rep," no-sum no-mean-imputation --bfile ",cur.dir,eth[i],"/all_chr.tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out ",cur.dir,eth[i],"/prs/prs_",k,"_rho_",l,"_size_",m))
    }
    
  }
  
  
}
time_prs1_end = proc.time()

#find the optimal r2 performance based on testing data
load(paste0(out.dir,"y_test.rdata"))
load(paste0(out.dir,"y_test.rdata"))
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
