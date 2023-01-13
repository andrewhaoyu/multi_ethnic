args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
#v = as.numeric(args[[3]])
#j = as.numeric(args[[3]])
#m = as.numeric(args[[3]])
#i_rep = as.numeric(args[[4]])
#i1 = as.numeric(args[[5]])

eth <- c("EUR","AFR","AMR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
trait = trait_vec[l] 

library(data.table)
library(tidyverse)

setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/GLGC_cleaned/"
kg.dir = "/data/zhangh24/KGref_MEGA/GRCh37/"
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')

#####
# Load library
library(XPASS)
library(data.table)
library(RhpcBLASctl)
#blas_set_num_threads(8)

file_out = paste0("/data/zhangh24/multi_ethnic/result/GLGC/prs/XPASS/",eth[i],"/",trait,"")
prs_mat = fread(paste0(file_out,"_PRS.sscore"))
colnames(prs_mat)[2] = "id"
pheno.dir = "/data/zhangh24/multi_ethnic/data/UKBB/phenotype/"
pheno_tuning = as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",eth[i],"_tuning.txt")))
pheno_tuning = pheno_tuning[,1:2]
covar <- as.data.frame(fread(paste0(pheno.dir,"/covariates/tuning+validation/",eth[i],"_all_data.txt")))
pheno_tuning <- left_join(pheno_tuning, covar)
colnames(pheno_tuning) = c('id','y','sex','age',paste0('pc',1:10))
pheno_tuning_com = pheno_tuning[complete.cases(pheno_tuning$y),]
pheno_tuning = left_join(pheno_tuning_com,prs_mat,by = "id")
pheno_vad = as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",eth[i],"_validation.txt")))
pheno_vad = pheno_vad[,1:2]
pheno_vad <- left_join(pheno_vad, covar)
colnames(pheno_vad) = c('id','y','sex','age',paste0('pc',1:10))
pheno_vad_com = pheno_vad[complete.cases(pheno_vad$y),]
pheno_vad = left_join(pheno_vad_com,prs_mat,by = "id")
#prs-csx doesn't require tuning
#use all data
pheno = rbind(pheno_tuning,pheno_vad)
model.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno)
y = model.null$residual
prs = pheno[,"SCORE1_SUM"]
model = lm(y~prs)
r2_xpass = summary(model)$r.square

data = data.frame(y = y, x = prs)
R2Boot = function(data,indices){
  boot_data = data[indices, ]
  model = lm(y ~ x, data = boot_data)
  result = summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 = boot(data = data, statistic = R2Boot, R = 10000)

ci_result = boot.ci(boot_r2, type = "bca")



r2_xpass = data.frame(eth = eth[i],
                       trait = trait_vec[l],
                       method = "XPASS",
                       r2 = r2,
                       r2_low = ci_result$bca[4],
                       r2_high = ci_result$bca[5]
)

out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/XPASS/",eth[i],"/",trait,"/")
save(r2_xpass, file = paste0(out.dir, "XPASS.result"))
# system(paste0("cp ",cur.dir,eth[i],"/all_chr_test.mega.* ",temp.dir))
# 
# # predict
# ref_gene_pred = paste0(temp.dir,"/all_chr_test.mega")
# # output file prefix
# load(paste0(file_out, "_param.RData"))
# mu = fit_bbj$mu
# # calculate PRS using plink (optional, comment the following out if not needed)
# mu = mu[, c("SNP", "A1", "mu1", "mu2", "mu_XPASS1", "mu_XPASS2")]
# write.table(mu,file = paste0(temp.dir,"prs_prep"),col.names = T,row.names = F,quote=F)
# res = system(paste0("/data/zhangh24/software/plink2_alpha ",
#                     "--score-col-nums 3,4,5,6 --threads 2 ",
#                     "--score ",temp.dir,"prs_prep cols=+scoresums,-scoreavgs header no-mean-imputation ",
#                     "--bfile ",ref_gene_pred,
#                     " --out ",file_out,"_PRS"))
# 




