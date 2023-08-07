#CT-SLEB for AOU data
#load the p-value results
args = commandArgs(trailingOnly = T)
#i represent ethnic group
#l represent trait

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
library(data.table)
library(dplyr)
library("CTSLEB", lib.loc = "/home/zhangh24/R/4.2/library")
eth_vec <- c("EUR","AFR","AMR")
trait_vec <-c("height","bmi")
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
data.dir = "/data/zhangh24/multi_ethnic/data/AOU_cleaned/"
out.dir = paste0("/data/zhangh24/multi_ethnic/result/AOU/clumping_result/PT/",eth,"/",trait,"/")

#load EUR data
sum_eur = as.data.frame(fread(paste0(data.dir,"EUR/",trait,"_update.txt"),header=T))


sum_eur = sum_eur %>% 
  select(rsID, CHR, pos37, BETA, SE, A1, P) %>% 
  rename(SNP = rsID, BP = pos37)
#load target population data
sum_tar = as.data.frame(fread(paste0(data.dir,eth,"/",trait,"_update.txt"),header=T))
sum_tar = sum_tar %>% 
  select(rsID, CHR, pos37, BETA, SE, A1, P) %>% 
  rename(SNP = rsID, BP = pos37)
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

load(paste0(out.dir, "all_SNP_set.rdata"))

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
eth_vec <- c("EUR","AFR","AMR")

eth_other = setdiff(eth_vec, eth_vec[i])

sum_other_list = list()
for(i_eth in 1:length(eth_other)){
  
  
  sum_temp = as.data.frame(fread(paste0(data.dir,eth_other[i_eth],"/",trait,"_update.txt"),header=T))
  sum_other_list[[i_eth]] = sum_temp %>% 
    select(rsID, CHR, pos37, BETA, SE, A1, P) %>% 
    rename(SNP = rsID, BP = pos37)
  
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


out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/AOU/prs/CT_SLEB_all/",eth,"/",trait,"/")


prs_list = list()
temp = 1
#take the column name of different clumping parameters
names = colnames(score_file)[3:ncol(score_file)]
for(k1 in 1:length(pthres)){
  for(k2 in 1:length(pthres)){
    #the --score file cols=+scoresums,-scoreavgs command in plink2 computes PRS as G*beta
    #If you compute PRS by chromosome, you need to sum the PRS scores for all chromosomes. 
    #load PRS for SNPs with p_target<p_thres[k2]|p_eur<p_thres[k1] 
    prs_temp = fread(paste0(out.dir.prs,"eb_prs_p_other_",k1,".p_tar_",k2,".sscore"))
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
prs_all_clean = prs_score %>% 
  select(-all_of(drop))

library(ranger)
library(SuperLearner)
#choose the prediction algorithms
SL.libray <- c(
  "SL.glmnet",
  "SL.ridge"
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
r2_ctsleb_true <- summary(model)$r.square






sl_fit = sl
#algorithm weight
alg_weights <- sl_fit$coef
#glmnet
glmnet_obj = sl_fit$fitLibrary$SL.glmnet$object
best_lambda <- glmnet_obj$lambda[which.min(glmnet_obj$cvm)]
glmnet_coefs <- coef(glmnet_obj, s = best_lambda)
#ridge
ridge_coefs = sl_fit$fitLibrary$SL.ridge_All$bestCoef
#final
final_coefs <- alg_weights[1] * glmnet_coefs + alg_weights[2] * ridge_coefs
#remove the intercept
final_coefs = final_coefs[2:nrow(final_coefs),]
#remove weight 0 coefficients
final_coefs = final_coefs[final_coefs!=0]
#clean ` in final_coefs
updated_name = gsub("`","",names(final_coefs))
names(final_coefs) = updated_name

modified_weights_names <- gsub("_p_other_.*", "", names(final_coefs))

# 1. Extract Relevant Columns from Score File
subset_score_columns <- score_file[, colnames(score_file) %in% c("SNP","A1",modified_weights_names), drop = FALSE]
p_value_file = score_file %>% 
  select(SNP) %>% 
  left_join(unique_infor_post) %>% 
  select(SNP,P,P_other)

# Function to extract p-value thresholds from weight names
extract_pvals <- function(weight_name) {
  p_other <- as.numeric(gsub(".*_p_other_([0-9e.-]+)_.*", "\\1", weight_name))
  p_tar <- as.numeric(gsub(".*_p_tar_([0-9e.-]+)$", "\\1", weight_name))
  return(list(p_other = p_other, p_tar = p_tar))
}


# Initialize an empty matrix to store results
results <- matrix(0, nrow = nrow(subset_score_columns), ncol = 1)
rownames(results) <- subset_score_columns$SNP
colnames(results) <- "final_weighted_score"

# Loop through each weight and apply the operations
for (weight_name in names(final_coefs)) {
  pvals <- extract_pvals(weight_name)
  score_colname <- gsub("_p_other_.*", "", weight_name)
  
  
  # Filter SNPs based on p-value criteria
  selected_snps <- na.omit(p_value_file$SNP[p_value_file$P <= pvals$p_tar | p_value_file$P_other <= pvals$p_other])
  
  # Multiply scores with the weight
  idx <- match(selected_snps, subset_score_columns$SNP)
  results[selected_snps, "final_weighted_score"] <- results[selected_snps, "final_weighted_score"] + 
    subset_score_columns[idx, score_colname] * final_coefs[weight_name]
}


#verify the pgs
prs_coef = cbind(score_file[,c("SNP","A1")],results)
write.table(prs_coef,file = paste0(temp.dir,"score_file_final_test"),row.names = F,col.names = F,quote=F)

res = system(paste0(soft.dir,"plink2_alpha ",
                    "--score-col-nums 3 ",
                    "--score ",temp.dir,"score_file_final_test cols=+scoresums,-scoreavgs ",
                    "--bfile ",temp.dir,"ukb/all_chr ",
                    "--out ",temp.dir,"eb_prs_final"))

prs_temp = fread(paste0(temp.dir,"eb_prs_final.sscore"))
# times (2*number of SNPs)
prs_mat = prs_temp[,c(1,2,5)]
colnames(prs_mat)[2] = "id"

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
prs_tun = pheno_tuning[,"SCORE1_SUM"]
pheno_vad = as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",eth,"_validation.txt")))
pheno_vad = pheno_vad[,1:2]
pheno_vad <- left_join(pheno_vad, covar)
colnames(pheno_vad) = c('id','y','sex','age',paste0('pc',1:10))
pheno_vad_com = pheno_vad[complete.cases(pheno_vad$y),]
pheno_vad = left_join(pheno_vad_com,prs_mat,by = "id")
model.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_vad)
y_vad = model.null$residual
prs_vad = pheno_vad[,"SCORE1_SUM"]
model <- lm(y_vad~prs_vad)
r2_ctsleb <- summary(model)$r.square
print(r2_ctsleb)
save(r2_ctsleb, file = paste0(out.dir, "CTSLEB_all_pgs.result"))



sum_tar = as.data.frame(fread(paste0(data.dir,eth,"/",trait,"_update.txt"),header=T))
sum_tar = sum_tar %>% 
  select(rsID, CHR, pos37, A2)

#prepare PRS for PGS catalog format
prs_select = left_join(prs_coef,sum_tar,by = c("SNP"="rsID")) %>% 
  rename(rsID = SNP,
         chr_name = CHR,
         chr_position = pos37,
         effect_allele = A1,
         other_allele = A2,
         effect_weight = final_weighted_score) %>% 
  select(rsID, chr_name, chr_position,
         effect_allele, other_allele, effect_weight)

out_filename = paste0("/data/zhangh24/multi_ethnic/result/AOU/pgs_catalog/",
                      trait,"_",eth,"_","CTSLEB.txt.gz")
write.table(prs_select, file = gzfile(out_filename), sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)
