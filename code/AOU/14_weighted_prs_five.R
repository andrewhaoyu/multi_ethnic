#LD_clumping for ARIC data
#load the p-value results
args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
#j = as.numeric(args[[3]])
#m = as.numeric(args[[3]])
#i_rep = as.numeric(args[[4]])
#i1 = as.numeric(args[[5]])
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)

eth_vec <- c("EUR","AFR","AMR")
trait_vec <-c("height","bmi")
eth = eth_vec[i]
trait = trait_vec[l]
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
# kg.dir = "/data/zhangh24/KGref_MEGA/GRCh37/"
# system(paste0("cp ",kg.dir,eth,"/all_chr.bed ",temp.dir,eth,"all_chr.bed"))
# system(paste0("cp ",kg.dir,eth,"/all_chr.bim ",temp.dir,eth,"all_chr.bim"))
# system(paste0("cp ",kg.dir,eth,"/all_chr.fam ",temp.dir,eth,"all_chr.fam"))

#load EUR results
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/AOU/"
#out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/EUR/",trait,"/")

load(paste0("/data/zhangh24/multi_ethnic/result/AOU/clumping_result/PT/EUR/",trait,"/CT.result"))
#find best index
idx <- which.max(ct.result[[2]])
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
data.dir = "/data/zhangh24/multi_ethnic/data/AOU_cleaned/"
#########prepare summary statistics#################
#prepare the coefficients for best CT PRS in other ans
for(i_eth in 1:length(eth_vec)){
  load(paste0("/data/zhangh24/multi_ethnic/result/AOU/clumping_result/PT/",eth_vec[i_eth],"/",trait,"/CT.result"))
  #find best index
  idx <- which.max(ct.result[[2]])
  pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
  
  #########prepare summary statistics#################
  #load eur summary statistics
  sum.data = as.data.frame(fread(paste0(data.dir,eth_vec[i_eth],"/",trait,"_update.txt"),header=T))
  sum.data.assoc = sum.data %>% 
    mutate(P = as.numeric(P)) %>% 
    rename(SNP = rsID,
           BP = pos37) %>% 
    select(CHR,SNP,BP,A1,BETA,P) 
  
  pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
  LD <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/AOU/clumping_result/PT/",
                                   eth_vec[i_eth],
                                   "/",trait,"/LD_clump.clumped")))
  clump.snp <- LD[,3,drop=F]
  #prepare the prs files with best EUR SNP and coefficients
  prs.all <- left_join(clump.snp,sum.data.assoc) %>% 
    filter(P <= pthres[idx])
  n_pthres <- length(pthres)
  
  #q_range = data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres),stringsAsFactors = F)
  prs.file = prs.all[,c("SNP","A1","BETA")]
  write.table(prs.file,file = paste0(temp.dir,"prs_coeff_",eth_vec[i_eth]),col.names = T,row.names = F,quote=F)
  
}

#########PRS calculation#################
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
out.dir = paste0("/data/zhangh24/multi_ethnic/result/AOU/weighted_prs/",eth,"/",trait,"/")
temp.dir = paste0('/lscratch/',sid,'/test/')
system(paste0("mkdir ",temp.dir,"ukb"))
geno.data = paste0("/data/zhangh24/multi_ethnic/data/UKBB/genotype/all_data/")
system(paste0("cp ",geno.data,eth,"/all_chr.bed ",temp.dir,"ukb/all_chr.bed"))
system(paste0("cp ",geno.data,eth,"/all_chr.bim ",temp.dir,"ukb/all_chr.bim"))
system(paste0("cp ",geno.data,eth,"/all_chr.fam ",temp.dir,"ukb/all_chr.fam"))
for(i_eth in 1:3){
  res <- system(paste0("/data/zhangh24/software/plink2_alpha ",
                       "--threads 2 ",
                       "--score ",temp.dir,"prs_coeff_",eth_vec[i_eth]," cols=+scoresums,-scoreavgs header no-mean-imputation ",
                       "--bfile ",temp.dir,"ukb/all_chr --out ",temp.dir,"prs_",eth_vec[i_eth]))
  
}

system(paste0("ls ",temp.dir,""))
#load prs
prs_temp = fread(paste0(temp.dir,"prs_EUR.sscore"))
prs_score = matrix(0, nrow(prs_temp), length(eth_vec))
for(i_eth in 1:length(eth_vec)){
  prs_temp = as.data.frame(fread(paste0(temp.dir,"prs_",eth_vec[i_eth],".sscore")))
  prs_score[, i_eth] = as.numeric(prs_temp[, "SCORE1_SUM"])
}
prs_infor = prs_temp[,1:4]
prs_mat = cbind(prs_infor, prs_score)
colnames(prs_mat)[5:7] = c(paste0("prs_", eth_vec))
colnames(prs_mat)[2] = "id"
#########R2 calculation################# 
pheno.dir = "/data/zhangh24/multi_ethnic/data/UKBB/phenotype/"
pheno_tuning = as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",eth,"_tuning.txt")))
pheno_tuning = pheno_tuning[,1:2]
covar <- as.data.frame(fread(paste0(pheno.dir,"/covariates/tuning+validation/",eth,"_all_data.txt")))
pheno_tuning <- left_join(pheno_tuning, covar)
colnames(pheno_tuning) = c('id','y','sex','age',paste0('pc',1:10))
pheno_tuning_com = pheno_tuning[complete.cases(pheno_tuning$y),]
pheno_tuning = left_join(pheno_tuning_com,prs_mat,by = "id")

pheno_vad = as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",eth,"_validation.txt")))
pheno_vad = pheno_vad[,1:2]
pheno_vad <- left_join(pheno_vad, covar)
colnames(pheno_vad) = c('id','y','sex','age',paste0('pc',1:10))
pheno_vad_com = pheno_vad[complete.cases(pheno_vad$y),]
pheno_vad = left_join(pheno_vad_com,prs_mat,by = "id")

#calculate R2 for each of the tuning dataset
prs_tun = as.matrix(pheno_tuning[,c(paste0("prs_", eth_vec))])
model.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_tuning)
model.prs <- lm(model.null$residual~prs_tun,data=pheno_tuning)
coeff = coefficients(model.prs)[2:4]

prs_vad = as.matrix(pheno_vad[,c(paste0("prs_", eth_vec))])
#evaluate on validation
model.vad.null  =  lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_vad)
prs  = prs_vad%*%coeff
model.vad.prs <- lm(model.vad.null$residual~prs,data=pheno_vad)
r2 = summary(model.vad.prs)$r.square

library(boot)

r2_function <- function(data, indices) {
  sample_data <- data[indices, ]
  model <- lm(y.vad ~ prs.vad, data = sample_data)
  return(summary(model)$r.square)
}
data <- data.frame(y.vad = model.vad.null$residual, prs.vad = prs)  # assuming y.vad and prs.vad are your data vectors
boot_r2 <- boot(data = data, statistic = r2_function, R = 10000)
ci_result = boot.ci(boot_r2, type = "perc")
ci_low= ci_result$percent[4]
ci_high = ci_result$percent[5]



result = data.frame(eth = eth,
                    trait = trait,
                    method = "Weighted PRS (C + T) five ancestries",
                    r2 = r2,
                    r2_low = ci_low,
                    r2_high = ci_high
)


save(result, file = paste0(out.dir, "weighted_prs_ct_three_ans.result"))

out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/AOU/prs/weightedprsall/",eth,"/",trait,"/")
prs_max_score = prs_score%*%coeff
prs_max = cbind(prs_infor, prs_max_score)
write.table(prs_max, file = paste0(out.dir.prs, "best_prs.sscore"),
            row.names = F,
            col.names = T,
            quote = F)



out_dir_boot = paste0("/data/zhangh24/multi_ethnic/result/AOU/boot_result/weighted_prs/",eth,"/",trait,"/")
boot_result = list(boot_r2,ci_result)
save(boot_result, file = paste0(out_dir_boot, "boot_result.rdata"))

