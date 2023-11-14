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

eth_vec <- c("EUR","AFR","AMR","EAS","SAS")
trait_vec <- c("HDL","LDL",
               "logTG",
               "TC")
eth = eth_vec[i]
trait = trait_vec[l]
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
kg.dir = "/data/zhangh24/KGref_MEGA/GRCh37/"
system(paste0("cp ",kg.dir,eth,"/all_chr.bed ",temp.dir,eth,"all_chr.bed"))
system(paste0("cp ",kg.dir,eth,"/all_chr.bim ",temp.dir,eth,"all_chr.bim"))
system(paste0("cp ",kg.dir,eth,"/all_chr.fam ",temp.dir,eth,"all_chr.fam"))

#load EUR results
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/GLGC_cleaned/"
#out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/EUR/",trait,"/")

load(paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/EUR/",trait,"/CT.result"))
#find best index
idx <- which.max(ct.result[[2]])
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
data.dir = "/data/zhangh24/multi_ethnic/data/GLGC_cleaned/"
#########prepare summary statistics#################
#load eur summary statistics
sum.data = as.data.frame(fread(paste0(data.dir,"EUR/",trait,".txt"),header=T))
sum.data.assoc = sum.data %>% 
  mutate(P = as.numeric(P)) %>% 
  rename(SNP = rsID,
         BP = POS_b37) %>% 
  select(CHR,SNP,BP,A1,BETA,P) 


pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
LD <- as.data.frame(fread(paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/EUR/",trait,"/LD_clump.clumped")))
clump.snp <- LD[,3,drop=F]
#prepare the prs files with best EUR SNP and coefficients
prs.all <- left_join(clump.snp,sum.data.assoc) %>% 
  filter(P <= pthres[idx])
n_pthres <- length(pthres)
#q_range = data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres),stringsAsFactors = F)
prs.file = prs.all[,c("SNP","A1","BETA")]
write.table(prs.file,file = paste0(temp.dir,"prs_coeff"),col.names = T,row.names = F,quote=F)

#########PRS calculation#################
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/BestEUR/",eth,"/",trait,"/")
temp.dir = paste0('/lscratch/',sid,'/test/')
system(paste0("mkdir ",temp.dir,"ukb"))
geno.data = paste0("/data/zhangh24/multi_ethnic/data/UKBB/genotype/all_data/")
system(paste0("cp ",geno.data,eth,"/all_chr.bed ",temp.dir,"ukb/all_chr.bed"))
system(paste0("cp ",geno.data,eth,"/all_chr.bim ",temp.dir,"ukb/all_chr.bim"))
system(paste0("cp ",geno.data,eth,"/all_chr.fam ",temp.dir,"ukb/all_chr.fam"))
res <- system(paste0("/data/zhangh24/software/plink2_alpha ",
                    "--threads 2 --score ",temp.dir,"prs_coeff cols=+scoresums,-scoreavgs header no-mean-imputation ",
                    "--bfile ",temp.dir,"ukb/all_chr --out ",temp.dir,"prs"))
system(paste0("ls ",temp.dir,""))
out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/GLGC/prs/EURPRS/",eth,"/",trait,"/")
system(paste0("cp ",temp.dir,"prs.sscore ",out.dir.prs))
#load prs
prs = fread(paste0(temp.dir,"prs.sscore"))
colnames(prs)[2] = "id"
#########R2 calculation################# 
pheno.dir = "/data/zhangh24/multi_ethnic/data/UKBB/phenotype/"
pheno = as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",eth,"_all_data.txt")))
pheno = pheno[,1:2]
covar <- as.data.frame(fread(paste0(pheno.dir,"/covariates/tuning+validation/",eth,"_all_data.txt")))
pheno <- left_join(pheno, covar)
colnames(pheno) = c('id','y','sex','age',paste0('pc',1:10))
pheno_com = pheno[complete.cases(pheno$y),]

#combine prs with pheno_all
pheno_all = left_join(pheno_com,prs,by = "id")

model.null  =  lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_all)
  prs = pheno_all[,paste0("SCORE1_SUM")]
  model.vad.prs <- lm(model.null$residual~prs,data=pheno_all)
  r2 = summary(model.vad.prs)$r.square
  
data = data.frame(y = model.null$residual, x = prs)
R2Boot = function(data,indices){
    boot_data = data[indices, ]
    model = lm(y ~ x, data = boot_data)
    result = summary(model)$r.square
    return(c(result))
  }
library(boot)
boot_r2 = boot(data = data, statistic = R2Boot, R = 20000)
  
ci_result = boot.ci(boot_r2, type = "bca")

r2.result = data.frame(eth = eth,
                       trait = trait,
                       method = "BestEUR_CT",
                       r2 = r2,
                       r2_low = ci_result$bca[4],
                       r2_high = ci_result$bca[5]
)

save(r2.result, file = paste0(out.dir, "BestEUR.result"))



out_dir_boot = out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/boot_result/BestEUR/",eth,"/",trait,"/")
boot_result = list(boot_r2,ci_result)
save(boot_result, file = paste0(out_dir_boot, "boot_result.rdata"))
