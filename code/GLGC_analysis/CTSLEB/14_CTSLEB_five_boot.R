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
setwd("/data/zhangh24/multi_ethnic/")

data.dir = "/data/zhangh24/multi_ethnic/data/GLGC_cleaned/"
out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/",eth,"/",trait,"/")
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)
out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/GLGC/prs/CT_SLEB_all/",eth,"/",trait,"/")
prs_list = list()
temp = 1
#take the column name of different clumping parameters
#names = colnames(score_file)[3:ncol(score_file)]
for(k1 in 1:length(pthres)){
  for(k2 in 1:length(pthres)){
    #the --score file cols=+scoresums,-scoreavgs command in plink2 computes PRS as G*beta
    #If you compute PRS by chromosome, you need to sum the PRS scores for all chromosomes. 
    #load PRS for SNPs with p_target<p_thres[k2]|p_eur<p_thres[k1] 
    prs_temp = fread(paste0(out.dir.prs,"eb_prs_p_other_",k1,".p_tar_",k2,".sscore"))
    # times (2*number of SNPs)
    prs_list[[temp]] = prs_temp[,5:ncol(prs_temp)]
    
    #colnames(prs_list[[temp]]) = paste0(names,"_","p_other_",pthres[k1],"_p_tar_",pthres[k2])
    temp = temp + 1
  }
}
prs_mat = as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
colnames(prs_mat)[2] = "id"
prs_score = prs_mat[,-c(1:2)]




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



data = data.frame(y = y_vad, x = y_pred[[1]])
R2Boot = function(data,indices){
  boot_data = data[indices, ]
  model = lm(y ~ x, data = boot_data)
  result = summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 = boot(data = data, statistic = R2Boot, R = 10000)

ci_result = boot.ci(boot_r2, type = "bca")

result = data.frame(eth = eth,
                    trait = trait,
                    method = "CT-SLEB (five ancestries)",
                    r2 = r2_ctsleb,
                    r2_low = ci_result$bca[4],
                    r2_high = ci_result$bca[5]
)



save(result, file = paste0(out.dir, "CTSLEB_all.result"))

#save the best prs
prs_max_score = predict(sl, prs_all_clean, onlySL = TRUE)[[1]]
prs_max = cbind(prs_temp[,1:4], prs_max_score)
write.table(prs_max, file = paste0(out.dir.prs, "best_prs.sscore"),
            row.names = F,
            col.names = T,
            quote = F)
#evaluate on validation


out_dir_boot = paste0("/data/zhangh24/multi_ethnic/result/GLGC/boot_result/CTSLEB/",eth,"/",trait,"/")
boot_result = list(boot_r2,ci_result)
save(boot_result, file = paste0(out_dir_boot, "boot_result.rdata"))


