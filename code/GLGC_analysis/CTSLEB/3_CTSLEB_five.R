#CT-SLEB for GLGC data
#load the p-value results
args = commandArgs(trailingOnly = T)
#i represent ethnic group
#l represent trait

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
library(data.table)
library(dplyr)
library(CTSLEB)
eth_vec <- c("EUR","AFR","AMR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
eth_EUR = "EUR"
eth = eth_vec[i]
trait = trait_vec[l]
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
kg.dir = "/data/zhangh24/KGref_MEGA/GRCh37/"
system(paste0("cp ",kg.dir,eth_EUR,"/all_chr.bed ",temp.dir,eth_EUR,"all_chr.bed"))
system(paste0("cp ",kg.dir,eth_EUR,"/all_chr.bim ",temp.dir,eth_EUR,"all_chr.bim"))
system(paste0("cp ",kg.dir,eth_EUR,"/all_chr.fam ",temp.dir,eth_EUR,"all_chr.fam"))

system(paste0("cp ",kg.dir,eth,"/all_chr.bed ",temp.dir,eth,"all_chr.bed"))
system(paste0("cp ",kg.dir,eth,"/all_chr.bim ",temp.dir,eth,"all_chr.bim"))
system(paste0("cp ",kg.dir,eth,"/all_chr.fam ",temp.dir,eth,"all_chr.fam"))
setwd("/data/zhangh24/multi_ethnic/")

data.dir = "/data/zhangh24/multi_ethnic/data/GLGC_cleaned/"
out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/",eth,"/",trait,"/")

#load EUR data
sum_eur = as.data.frame(fread(paste0(data.dir,"EUR/",trait,".txt"),header=T))


sum_eur = sum_eur %>% 
  select(rsID, CHR, POS_b37, BETA, SE, A1, P) %>% 
  rename(SNP = rsID, BP = POS_b37)
#load target population data
sum_tar = as.data.frame(fread(paste0(data.dir,eth,"/",trait,".txt"),header=T))
sum_tar = sum_tar %>% 
  select(rsID, CHR, POS_b37, BETA, SE, A1, P) %>% 
  rename(SNP = rsID, BP = POS_b37)
#align allels
sum_com <- AlignSum(sum_tar = sum_tar,
                    sum_other = sum_eur)

#split the SNPs into two groups
sum_com_split <- SplitSum(sum_com)
#sum_com_split is a list with two data frame
#the first data frame contains SNPs with p_eur < p_target, the p-value column is from p_eur.
sum_other_ref = sum_com_split[[1]]
#the second data frame contains target population-specific SNPs or p_eur < p_target. The p-value column is from p_target. 
sum_tar_ref = sum_com_split[[2]]




###################CLUMPING STEP started##############

# we use plink1.9 for the clumping purpose
# specify vector for clumping r square and base window size
#the clumping_window_ize = base_window_size/clumping_r_square so that lower clumping r2 can have larger clumping window size
r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)

write.table(sum_other_ref,paste0(temp.dir,"sum_other_ref"),col.names = T,row.names = F,quote=F)
write.table(sum_tar_ref,paste0(temp.dir,"sum_tar_ref"),col.names = T,row.names = F,quote=F)
# --clump-p1 determines the upper bound of p-value to be kept in the clumping. We set it as 1.
# --clump-r2 is the clumping r square
# --clump-kb is the clumping window size
snp_list = list()
temp = 1
soft.dir = "/data/zhangh24/software/"
snp_list = list()
temp = 1
for(r_ind in 1:length(r2_vec)){
  #create the window size given the clumping r2
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    pthr = 1.0
    r2thr = r2_vec[r_ind]
    kbpthr = wc_vec[w_ind]
    #for the first group, we perform clumping using EUR popultion as the reference   
    system(paste0(soft.dir,"plink2 ",
                  "--bfile ",temp.dir,eth_EUR,"all_chr ",
                  "--clump ",temp.dir,"sum_other_ref ",
                  "--clump-p1 ",pthr," ",
                  "--clump-r2 ",r2thr," ",
                  "--clump-kb ",kbpthr," ",
                  "--out ", temp.dir,"EUR_ref_CT_rind_",r_ind,"_wcind_",w_ind))
    #for the second group, we perform clumping using AFR population as the reference
    system(paste0(soft.dir,"plink2 ",
                  "--bfile ",temp.dir,eth,"all_chr ",
                  "--clump ",temp.dir,"sum_tar_ref ",
                  "--clump-p1 ",pthr," ",
                  "--clump-r2 ",r2thr," ",
                  "--clump-kb ",kbpthr," ",
                  "--out ", temp.dir,eth,"_ref_CT_rind_",r_ind,"_wcind_",w_ind))
    
    #combine the SNPs from the two clumping groups
    LD_EUR= fread(paste0(temp.dir,"EUR_ref_CT_rind_",r_ind,"_wcind_",w_ind,".clumped"))[,3,drop=F]
    LD_tar = fread(paste0(temp.dir,eth,"_ref_CT_rind_",r_ind,"_wcind_",w_ind,".clumped"))[,3,drop=F]
    LD  = rbind(LD_EUR,LD_tar)
    snp_list[[temp]] = LD
    names(snp_list[[temp]]) = paste0("clump_r2_",r2thr,"_ws_",kbpthr)
    temp = temp + 1
  }
}

################CLUMPING STEP FINISHED###############

############################PRS step###############
system(paste0("mkdir ",temp.dir,"ukb"))
geno.data = paste0("/data/zhangh24/multi_ethnic/data/UKBB/genotype/all_data/")
system(paste0("cp ",geno.data,eth,"/all_chr.bed ",temp.dir,"ukb/all_chr.bed"))
system(paste0("cp ",geno.data,eth,"/all_chr.bim ",temp.dir,"ukb/all_chr.bim"))
system(paste0("cp ",geno.data,eth,"/all_chr.fam ",temp.dir,"ukb/all_chr.fam"))

plink_file = PreparePlinkFile(snp_list,sum_com)

score_file = plink_file[[1]]
write.table(score_file,file = paste0(temp.dir,"score_file"),row.names = F,col.names = F,quote=F)
p_value_file = plink_file[[2]]
unique_infor = plink_file[[3]]
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
#create q-range file
q_range = CreateQRange(pthres)
head(q_range)
write.table(q_range,file = paste0(temp.dir,"q_range_file"),row.names = F,col.names = F,quote=F)

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
  #AFR_test_chr22 contains 20,000 subjects
  #we use the first 10,000 subjects as tuning dataset
  #we use the second 10,000 subjects as validation dataset      
  res = system(paste0(soft.dir,"plink2_alpha ",
                      "--q-score-range ",temp.dir,"q_range_file ",temp.dir,"p_value_file ",
                      "--score-col-nums 3-",n_col," ",
                      "--score ",temp.dir,"score_file cols=+scoresums,-scoreavgs ",
                      "--bfile ",temp.dir,"ukb/all_chr ",
                      "--out ",temp.dir,"prs_p_other_",k1))
  
}
prs_list = list()
temp = 1
#take the column name of different clumping parameters
names = colnames(score_file)[3:ncol(score_file)]
for(k1 in 1:length(pthres)){
  for(k2 in 1:length(pthres)){
    #the --score file cols=+scoresums,-scoreavgs command in plink2 computes PRS as G*beta
    #If you compute PRS by chromosome, you need to sum the PRS scores for all chromosomes. 
    #load PRS for SNPs with p_target<p_thres[k2]|p_eur<p_thres[k1] 
    prs_temp = fread(paste0(temp.dir,"prs_p_other_",k1,".p_tar_",k2,".sscore"))
    # times (2*number of SNPs)
    prs_list[[temp]] = prs_temp[,5:ncol(prs_temp)]
    
    colnames(prs_list[[temp]]) = paste0(names,"_","p_other_",pthres[k1],"_p_tar_",pthres[k2])
    temp = temp + 1
  }
}
prs_mat = as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
colnames(prs_mat)[2] = "id"



#find the best R-square among the all the PRSs to find candidate set
#we use this candidate for estimating covariance matrix for the prior distribution
#create prediction r2 vector to store r2 for different prs
n.total.prs = length(pthres)^2*length(r2_vec)*length(wc_base_vec)
prs_r2_vec_test = rep(0,n.total.prs)
#load the phenotype data for the tuning set
pheno.dir = "/data/zhangh24/multi_ethnic/data/UKBB/phenotype/"
pheno_tuning = as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",eth,"_tuning.txt")))
pheno_tuning = pheno_tuning[,1:2]
covar <- as.data.frame(fread(paste0(pheno.dir,"/covariates/tuning+validation/",eth,"_all_data.txt")))
pheno_tuning <- left_join(pheno_tuning, covar)
colnames(pheno_tuning) = c('id','y','sex','age',paste0('pc',1:10))
pheno_tuning_com = pheno_tuning[complete.cases(pheno_tuning$y),]
pheno_tuning = left_join(pheno_tuning_com,prs_mat,by = "id")
model.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_tuning)
y_tun = model.null$residual
prs_tun = pheno_tuning[, colnames(prs_mat)]
for(p_ind in 1:n.total.prs){
  #the first two columns of prs_tun are family id and individual id
  #prs starts from the third column
  model = lm(y_tun~prs_tun[,(2+p_ind)])
  prs_r2_vec_test[p_ind] = summary(model)$r.square
}
#+2 is due to the first two columns are family id and individual id
max_ind = which.max(prs_r2_vec_test)
print(colnames(prs_mat)[max_ind+2])

######################PRS step finished#############




#############EB step start############################

#Get the SNP set with the best performance in the CT step
snp_set_ind = colnames(prs_mat)[max_ind+2]
SNP_set = GetSNPSet(snp_set_ind,
                    score_file,
                    unique_infor)
save(SNP_set, file = paste0(out.dir, "SNP_set.rdata"))

##########five ancestries analyses###############
#########EB step###################
eth_vec <- c("EUR","AFR","AMR","EAS", "SAS")

eth_other = setdiff(eth_vec, eth_vec[i])

sum_other_list = list()
for(i_eth in 1:length(eth_other)){
  
  
  sum_temp = as.data.frame(fread(paste0(data.dir,eth_other[i_eth],"/",trait,".txt"),header=T))
  sum_other_list[[i_eth]] = sum_temp %>% 
    select(rsID, CHR, POS_b37, BETA, SE, A1, P) %>% 
    rename(SNP = rsID, BP = POS_b37)
  
}
other_ans_names = eth_other
sum_com <- AlignSumMulti(sum_tar = sum_tar,
                         sum_other_list = sum_other_list,
                         other_ans_names = other_ans_names)



unique_infor_post = EBpostMulti(unique_infor,SNP_set,
                                sum_com,other_ans_names)

eb_post_col_names = c("BETA_EB_target",paste0("BETA_EB_",other_ans_names[1]))
post_beta_mat = unique_infor_post %>% 
  select(all_of(eb_post_col_names))


plink_file_eb = PreparePlinkFileEB(snp_list,
                                   unique_infor_post,
                                   post_beta_mat)
score_file = plink_file_eb[[1]]
write.table(score_file,file = paste0(temp.dir,"score_file_eb"),row.names = F,col.names = F,quote=F)
#p_value_file description
#the second column is the p-values of SNPs from the GWAS of the target population
p_value_file = plink_file_eb[[2]]


p_value_file_temp = p_value_file
for(k1 in 1:length(pthres)){
  #keep al the SNPs with P_EUR less than pthres[k1] in the analyses
  idx <- which(unique_infor$P_other<=pthres[k1])
  p_value_file_temp$P[idx] = 0
  write.table(p_value_file_temp,file = paste0(temp.dir,"p_value_file"),col.names = F,row.names = F,quote=F)
  n_col = ncol(score_file)
  
  res = system(paste0(soft.dir,"plink2_alpha ",
                      "--q-score-range ",temp.dir,"q_range_file ",temp.dir,"p_value_file ",
                      "--score-col-nums 3-",n_col," ",
                      "--score ",temp.dir,"score_file_eb cols=+scoresums,-scoreavgs ",
                      "--bfile ",temp.dir,"ukb/all_chr ",
                      "--out ",temp.dir,"eb_prs_p_other_",k1))
  #the output of plink2 create 9 different files named as prs_p_other_k1.p_tar_k2.sscore
  #this output file contains 16 columns
  #the column 1-4 are: family ID, individual ID, 2*total number of SNPs in the PRS, the sum of allele count
  #column 5-16 are the PRS scores with SNP of p_target<p_thres[k2]|p_eur<p_thres[k1] for different combinations of r2-cutoff and base_window_size
}
#find best cutoff for EUR by using all data as tuning

prs_list = list()
temp = 1
#take the column name of different clumping parameters
names = colnames(score_file)[3:ncol(score_file)]
for(k1 in 1:length(pthres)){
  for(k2 in 1:length(pthres)){
    #the --score file cols=+scoresums,-scoreavgs command in plink2 computes PRS as G*beta
    #If you compute PRS by chromosome, you need to sum the PRS scores for all chromosomes. 
    #load PRS for SNPs with p_target<p_thres[k2]|p_eur<p_thres[k1] 
    prs_temp = fread(paste0(temp.dir,"eb_prs_p_other_",k1,".p_tar_",k2,".sscore"))
    # times (2*number of SNPs)
    prs_list[[temp]] = prs_temp[,5:ncol(prs_temp)]
    
    colnames(prs_list[[temp]]) = paste0(names,"_","p_other_",pthres[k1],"_p_tar_",pthres[k2])
    temp = temp + 1
  }
}
prs_mat = as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
colnames(prs_mat)[2] = "id"
prs_score = prs_mat[,-c(1:2)]
#move the prs to the prs file
out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/GLGC/prs/CT_SLEB_all/",eth,"/",trait,"/")
system(paste0("cp ",temp.dir,"eb_prs_p_other_*.sscore ",out.dir.prs))

#############EB step finish############################

#system(past)





############SL step start#############################
pheno.dir = "/data/zhangh24/multi_ethnic/data/UKBB/phenotype/"
pheno_tuning = as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",eth,"_tuning.txt")))
pheno_tuning = pheno_tuning[,1:2]
covar <- as.data.frame(fread(paste0(pheno.dir,"/covariates/tuning+validation/",eth,"_all_data.txt")))
pheno_tuning <- left_join(pheno_tuning, covar)
colnames(pheno_tuning) = c('id','y','sex','age',paste0('pc',1:10))
pheno_tuning_com = pheno_tuning[complete.cases(pheno_tuning$y),]
pheno_tuning = left_join(pheno_tuning_com,prs_mat,by = "id")
model.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_tuning)
y_tun = model.null$residual
prs_tun = pheno_tuning[,colnames(prs_score)]
pheno_vad = as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",eth,"_validation.txt")))
pheno_vad = pheno_vad[,1:2]
pheno_vad <- left_join(pheno_vad, covar)
colnames(pheno_vad) = c('id','y','sex','age',paste0('pc',1:10))
pheno_vad_com = pheno_vad[complete.cases(pheno_vad$y),]
pheno_vad = left_join(pheno_vad_com,prs_mat,by = "id")
model.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_vad)
y_vad = model.null$residual
prs_vad = pheno_vad[,colnames(prs_score)]
library(caret)

mtx = cor(prs_tun)
drop = findCorrelation(mtx,cutoff=0.98)
drop = names(prs_tun)[drop]

prs_tun_clean = prs_tun %>% 
  select(-all_of(drop))
prs_vad_clean = prs_vad %>% 
  select(-all_of(drop))
prs_all_clean = prs_score %>% 
  select(-all_of(drop))
library(ranger)
library(SuperLearner)
#choose the prediction algorithms
SL.libray <- c(
  "SL.glmnet",
  "SL.ridge",
  "SL.nnet"
  #"SL.bayesglm"
  #"SL.stepAIC"
  #"SL.xgboost"
  #"SL.randomForest"
  #"SL.ksvm",
  #"SL.bartMachine", 
  #"SL.kernelKnn",
  #"SL.rpartPrune", 
  #"SL.lm"
  #"SL.mean"
)
sl = SuperLearner(Y = y_tun, X = prs_tun_clean, family = gaussian(),
                  # For a real analysis we would use V = 10.
                  # V = 3,
                  SL.library = SL.libray)

y_pred <- predict(sl, prs_vad_clean, onlySL = TRUE)
#evaluate the CT-SLEB prs performance on the validation
model <- lm(y_vad~y_pred[[1]])
r2_ctsleb <- summary(model)$r.square





save(r2_ctsleb, file = paste0(out.dir, "CTSLEB_all.result"))

#save the best prs
prs_max_score = predict(sl, prs_all_clean, onlySL = TRUE)[[1]]
prs_max = cbind(prs_temp[1:4], prs_max_score)
write.table(prs_max, file = paste0(out.dir.prs, "best_prs.sscore"),
            row.names = F,
            col.names = T,
            quote = F)
#evaluate on validation


