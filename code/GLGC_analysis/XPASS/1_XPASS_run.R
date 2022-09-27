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
system(paste0("cp ",kg.dir,eth[1],"/all_chr.bed ",temp.dir,eth[1],"all_chr.bed"))
system(paste0("cp ",kg.dir,eth[1],"/all_chr.bim ",temp.dir,eth[1],"all_chr.bim"))
system(paste0("cp ",kg.dir,eth[1],"/all_chr.fam ",temp.dir,eth[1],"all_chr.fam"))
system(paste0("cp ",kg.dir,eth[1],"/all_chr_pca.eigenvec ",temp.dir,eth[1],"all_chr_pca.eigenvec"))
system(paste0("cp ",kg.dir,eth[i],"/all_chr.bed ",temp.dir,eth[i],"all_chr.bed"))
system(paste0("cp ",kg.dir,eth[i],"/all_chr.bim ",temp.dir,eth[i],"all_chr.bim"))
system(paste0("cp ",kg.dir,eth[i],"/all_chr.fam ",temp.dir,eth[i],"all_chr.fam"))
system(paste0("cp ",kg.dir,eth[i],"/all_chr_pca.eigenvec ",temp.dir,eth[i],"all_chr_pca.eigenvec"))
#####
# Load library
library(XPASS)
library(data.table)
library(RhpcBLASctl)
#blas_set_num_threads(8)

#####
# parse the command line argument
#eval(parse(text=paste(commandArgs(trailingOnly = TRUE), collapse=";")))

### if not using command line, then please specify the parameters below.
# race = "AFR" # target population
# size = 1 # size for target population
# GA = 1
# rho = 1
# rep = 1
#####
# input files
summary_EUR = as.data.frame(fread(paste0(data.dir,"EUR/",trait,".txt"),header=T))
summary_target = as.data.frame(fread(paste0(data.dir,eth[i],"/",trait,".txt"),header=T))



#prepare the EUR data for XPASS format
summary_EUR_XPASS = summary_EUR %>% 
  mutate(Z = BETA/SE,
         SNP = rsID) %>% 
  select(SNP, N, Z, A1, A2) 
write.table(summary_EUR_XPASS, file = paste0(temp.dir,eth[1],"summary_eur"),
            row.names = F, col.names = T, quote = T)
#summary_EUR = paste0("/dcs04/nilanjan/data/wlu/XPASS/data/sumdata/EUR/sumdata-rho",rho,'-size4-rep',rep,'-GA',GA,'.txt') # auxilliary

summary_target_XPASS = summary_target %>% 
  mutate(Z = BETA/SE,
         SNP = rsID) %>% 
  select(SNP, N, Z, A1, A2) 
write.table(summary_target_XPASS, file = paste0(temp.dir,eth[i],"summary_tar"),
            row.names = F, col.names = T, quote = T)
path_to_summary_eur = paste0(temp.dir,eth[1],"summary_eur")
path_to_summary_tar = paste0(temp.dir,eth[i],"summary_tar")
#summary_target = paste0("/dcs04/nilanjan/data/wlu/XPASS/data/sumdata/",race,'/sumdata-rho',rho,'-size',size,'-rep',rep,'-GA',GA,'.txt') # target
# auxilliary
ref_gene_EUR = paste0(temp.dir,eth[1],"all_chr")
# target
ref_gene_target = paste0(temp.dir,eth[i],"all_chr")
file_out = paste0("/data/zhangh24/multi_ethnic/result/GLGC/prs/XPASS/",eth[i],"/",trait,"")

eur_ref_cov = as.data.frame(fread(paste0(temp.dir,eth[1],"all_chr_pca.eigenvec")))
eur_ref_cov = eur_ref_cov[,3:22]
tar_ref_cov = as.data.frame(fread(paste0(temp.dir,eth[i],"all_chr_pca.eigenvec")))
tar_ref_cov = tar_ref_cov[,3:22]
write.table(eur_ref_cov, file = paste0(temp.dir,eth[1],"eur_ref_cov"),
            row.names = F, col.names = F, quote = F)

write.table(tar_ref_cov, file = paste0(temp.dir,eth[i],"tar_ref_cov"),
            row.names = F, col.names = F, quote = F)

path_to_cov_EUR = paste0(temp.dir,eth[1],"eur_ref_cov")
path_to_cov_tar = paste0(temp.dir,eth[i],"tar_ref_cov")

#####
# XPASS
fit_bbj <- XPASS(file_z1 = path_to_summary_tar, file_z2 = path_to_summary_eur,
                 file_ref1 = ref_gene_target, file_ref2 = ref_gene_EUR,
                 file_cov1 = path_to_cov_tar,file_cov2 = path_to_cov_EUR,
                 # file_predGeno = ref_gene_pred, compPRS = T,
                 pop = "EUR", sd_method="LD_block", compPosMean = T,
                 file_out = file_out)
save(fit_bbj, file=paste0(file_out, "_param.RData"))
print("XPASS finished run")

system(paste0("mkdir ",temp.dir,"ukb"))
geno.data = paste0("/data/zhangh24/multi_ethnic/data/UKBB/genotype/all_data/")
system(paste0("cp ",geno.data,eth[i],"/all_chr.bed ",temp.dir,"ukb/all_chr.bed"))
system(paste0("cp ",geno.data,eth[i],"/all_chr.bim ",temp.dir,"ukb/all_chr.bim"))
system(paste0("cp ",geno.data,eth[i],"/all_chr.fam ",temp.dir,"ukb/all_chr.fam"))



mu = fit_bbj$mu
mu = mu[, c("SNP", "A1", "mu_XPASS1")]
write.table(mu,file = paste0(temp.dir,"prs_prep"),col.names = T,row.names = F,quote=F)

ref_gene_pred = paste0(temp.dir,"ukb/all_chr")
res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                    "--score-col-nums 3 --threads 2 ",
                    "--score ",temp.dir,"prs_prep cols=+scoresums,-scoreavgs header no-mean-imputation ",
                    "--bfile ",ref_gene_pred,
                    " --out ",file_out,"_PRS"))

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




