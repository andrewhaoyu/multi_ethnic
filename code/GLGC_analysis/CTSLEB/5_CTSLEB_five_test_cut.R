#CT-SLEB for GLGC data
#load the p-value results
args = commandArgs(trailingOnly = T)
#i represent ethnic group
#l represent trait

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
z_ind = as.numeric(args[[3]])
library(data.table)
library(dplyr)
library(CTSLEB)
pthres <- c(Inf,1E-10,5E-08,5E-05,1.0)
#when pthres is Inf; all SNPs will be included in the EB step
#when ptres is 1; no SNP will be included in EB step; just use target populatoin coeff
z_cut = -qnorm(pthres[z_ind]/2)
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
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
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

GetSNPSet = function(snp_set_ind,
                     score_file,
                     unique_infor){
  str_temp = strsplit(snp_set_ind,"_")
  r2 = str_temp[[1]][3]
  ws = str_temp[[1]][5]
  p_other_cutoff = as.numeric(str_temp[[1]][[8]])
  p_tar_cutoff = as.numeric(str_temp[[1]][[11]])
  
  score_file_name = paste0("clump_r2_",r2,"_ws_",ws)
  idx = which(colnames(score_file)==score_file_name)
  
  #take the post-clupmping SNPs with particular r2-cutoff and window size by removing SNPs with coefficients 0
  LD = score_file[score_file[,idx]!=0,"SNP",drop=F]
  #take the p_value cutoff
  LD_infor = left_join(LD,unique_infor,by="SNP")
  #take the SNP with p_eur < p_other_cutoff | p_tar < p_tar_cutoff
  snp_ind = which(LD_infor$P_other<=p_other_cutoff|
                    LD_infor$P<=p_tar_cutoff)
  SNP = LD_infor[snp_ind,]
  return(SNP)
}



#############EB step start############################

#Get the SNP set with the best performance in the CT step
snp_set_ind = colnames(prs_mat)[max_ind+2]
SNP_set = GetSNPSet(snp_set_ind,
                    score_file,
                    unique_infor)

#save(SNP_set, file = paste0(out.dir, "SNP_set.rdata"))

AlignSumMulti = function(sum_tar,sum_other_list,
                         other_ans_names){
  coeff_other_list = list()
  for(i in 1:length(other_ans_names)){
    sum_other_temp  = sum_other_list[[i]]
    sum_com_temp <- AlignSum(sum_tar = sum_tar,
                             sum_other = sum_other_temp)
    #a temporary matrix to save the aligned cofficients for the target population
    coeff_other = sum_com_temp[,c("BETA_other","SE_other","P_other")]
    colnames(coeff_other) = paste0(c("BETA_","SE_","P_"),other_ans_names[i])
    coeff_other_list[[i]] = coeff_other
  }
  coeff_mat = bind_cols(coeff_other_list)
  
  sum_com = cbind(sum_tar,coeff_mat)
  return(sum_com)
}

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

EstimatePriorMultiUpdate <- function(SNP_set,other_ans_names,
                               sum_com, z_cut = z_cut){
  if(is.nan(z_cut)){
    z_cut = -qnorm(5E-08/2) #keep geno-wide significant snps out of analyse
    }
  SNP_set_select = SNP_set %>%
    select(SNP)
  #align the summary statistics with the SNP_set from CT
  SNP_set_align = left_join(SNP_set_select,sum_com,
                            by="SNP")
  n_ans = length(other_ans_names)
  
  #create column names to select the BETA_ and SE_ columns
  col_names_beta = c("BETA",paste0("BETA_",other_ans_names))
  col_names_se = c("SE",paste0("SE_",other_ans_names))
  beta_mat = SNP_set_align %>%
    select(all_of(col_names_beta))
  se_mat = SNP_set_align %>%
    select(all_of(col_names_se))
  #\hat_{u}_kl|u_kl ~ N(u_kl,1/N_l)
  #where u_kl is the underlying effect size for the kth SNP of lth population
  #N_l is the sample size
  #since the Bayesian algorithm is applied on the standardized effect-size scale
  #it's equivalent to applying the Bayes rule on z-statistics scale
  #the advantage of z-statistics scale is the covariance matrix is identity
  #it can make the computation faster.
  z_mat = beta_mat/se_mat
 
  eb_idx = which(apply(z_mat, 1, function(r) any(abs(r)>=z_cut, na.rm = T))==F)
  if(length(eb_idx)>0){
    #only use SNPs without large effects to estimate prior
    z_mat_clean = z_mat[eb_idx,]
    
    z_mat_clean <-na.omit(z_mat_clean)
    p = ncol(z_mat_clean)
    
    prior_mat <- cov(z_mat_clean, use='pairwise')-diag(p)
    colnames(prior_mat) = c("Z_tar",paste0("Z_",other_ans_names))
    
  }else{
    prior_mat = NULL
  }
  return(prior_mat)
}

EBpostMultiUpdate <- function(unique_infor,SNP_set,
                        sum_com,other_ans_names,z_cut){
  
  #prior_sigma = prior_result[[1]]
  #eb_idx = prior_result[[2]]
  SNP_set_select = unique_infor %>%
    select(SNP)
  #align the summary statistics with the SNP_set from CT
  SNP_set_align = left_join(SNP_set_select,sum_com,
                            by="SNP")
  
  
  #create column names to select the BETA_ and SE_ columns
  col_names_beta = c("BETA",paste0("BETA_",other_ans_names))
  col_names_se = c("SE",paste0("SE_",other_ans_names))
  beta_mat = SNP_set_align %>%
    select(all_of(col_names_beta))
  se_mat = SNP_set_align %>%
    select(all_of(col_names_se))
  
  z_mat = as.matrix(beta_mat/se_mat)
  
  prior_sigma = EstimatePriorMultiUpdate(SNP_set,other_ans_names,
                                          sum_com,z_cut) 
  if(is.null(prior_sigma)){

    colnames(beta_mat) = c("BETA_EB_target",paste0("BETA_EB_",other_ans_names))
    eb_beta_names = colnames(beta_mat)
    unique_infor_EB =cbind(unique_infor,beta_mat) %>%
      select(SNP,A1,all_of(eb_beta_names),P,P_other)
    return(unique_infor_EB)
  }else{

    #filter out the SNPs with large effects out of EB step
    eb_idx = which(apply(z_mat, 1, function(r) any(abs(r)>=z_cut, na.rm = T))==F)
    non_eb_idx = which(apply(z_mat, 1, function(r) any(abs(r)>=z_cut, na.rm = T))==T)
    z_mat_post = as.matrix(z_mat)
    col_names_beta = c("Z",paste0("Z_",other_ans_names))
    p <- ncol(z_mat)
    
    post_sigma = solve(solve(prior_sigma)+diag(p))
    
    #if you don't want to use EB procedure, you can set z_cut to be 0, then eb_idx will be NULL
    #under this specicial scenaior, we will just use the original beta from the target population
    if(length(eb_idx)==0){
      colnames(beta_mat) = c("BETA_EB_target",paste0("BETA_EB_",other_ans_names))
      eb_beta_names = colnames(beta_mat)
      unique_infor_EB =cbind(unique_infor,beta_mat) %>%
        select(SNP,A1,all_of(eb_beta_names),P,P_other)
    }else{
      for(k in 1:length(eb_idx)){
        if(k%%10000==0){print(paste0(k," SNPs completed"))}
        z_temp =z_mat[eb_idx[k],]
        
        #find out nonmissing component
        
        idx <- which(!is.na(z_temp))
        if(length(idx)<p){
          z_temp <- z_temp[idx]
          
          post_sigma_temp = post_sigma[idx,idx,drop=F]
          z_post = post_sigma_temp%*%z_temp
        }else{
          z_post =post_sigma%*%z_temp
        }
        
        z_mat_post[eb_idx[k],idx] = z_post
      }
      beta_mat_post = z_mat_post*se_mat
      colnames(beta_mat_post) = c("BETA_EB_target",paste0("BETA_EB_",other_ans_names))
      eb_beta_names = colnames(beta_mat_post)
      unique_infor_EB =cbind(unique_infor,beta_mat_post) %>%
        select(SNP,A1,all_of(eb_beta_names),P,P_other)
      
    }
    return(unique_infor_EB)
  }
 
  
}



FindUniqueSNP = function(snp,
                         sum_com){
  unique_id = unique(rbindlist(snp_list,use.name =FALSE))
  names(unique_id) = "SNP"
  #align the regression coefficients for these SNPs from the sum stat
  unique_infor = left_join(unique_id,sum_com,by="SNP")
  return(unique_infor)
}
unique_infor = FindUniqueSNP(snp_list,sum_com)

unique_infor_post = EBpostMultiUpdate(unique_infor,SNP_set,
                                sum_com,other_ans_names,z_cut)

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

#############EB step finish############################







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





save(r2_ctsleb, file = paste0(out.dir, "CTSLEB_all_test_",z_ind,".result"))


