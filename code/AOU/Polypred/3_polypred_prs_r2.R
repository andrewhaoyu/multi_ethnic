args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])

# for(i in 2:3){
#   for(l in 1:2){
library(dplyr)
library(data.table)
eth <- c("EUR","AFR","AMR")
trait <- c("height","bmi")

out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/AOU/prs/PRSCSX/",eth[i],"/",trait,"/")
phi = c("1e+00","1e-02","1e-04","1e-06")

sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')

system(paste0("mkdir ",temp.dir,"ukb"))
geno.data = paste0("/data/zhangh24/multi_ethnic/data/UKBB/genotype/all_data/")
system(paste0("cp ",geno.data,eth[i],"/all_chr.bed ",temp.dir,"ukb/all_chr.bed"))
system(paste0("cp ",geno.data,eth[i],"/all_chr.bim ",temp.dir,"ukb/all_chr.bim"))
system(paste0("cp ",geno.data,eth[i],"/all_chr.fam ",temp.dir,"ukb/all_chr.fam"))
ref_gene_pred = paste0(temp.dir,"ukb/all_chr")
#calculate the PRS for SBayesR based on EUR 
file_dir = paste0("/data/zhangh24/multi_ethnic/result/AOU/prs/polypred/",eth[1],"/",trait[l])
#coefs are based on EUR population
files = dir(file_dir, pattern = ".snpRes", full.names = T)
eur_coef = paste0("/data/zhangh24/multi_ethnic/result/AOU/prs/polypred/",eth[1],"/",trait[l],"/SBayesR.snpRes")

#output file needs to be the target population


if(eur_coef%in%files==T){
  #load result
  sbayes_result = fread(paste0(eur_coef),header = T)
  #if the job converges, then nrow(file_coef) > 0
  if(nrow(sbayes_result)>0){
      
    
    assoc = sbayes_result %>% 
      rename(BETA = A1Effect,
             SNP = Name) %>% 
      select(SNP, A1, BETA)
    write.table(assoc,file = paste0(temp.dir,"prs_prep"),col.names = T,row.names = F,quote=F)
    #calculate the PRS for SBayesR 
    
    res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                        "--score-col-nums 3 --threads 2 ",
                        "--score ",temp.dir,"prs_prep cols=+scoresums,-scoreavgs header no-mean-imputation ",
                        "--bfile ",ref_gene_pred,
                        " --out ",temp.dir,"PRS_SBayesR_EUR"))
  }
  
}



#calculate the PRS for SBayesR based on target population


file_dir = paste0("/data/zhangh24/multi_ethnic/result/AOU/prs/polypred/",eth[i],"/",trait[l])
#coefs are based on EUR population
files = dir(file_dir, pattern = ".snpRes", full.names = T)
tar_coef = paste0("/data/zhangh24/multi_ethnic/result/AOU/prs/polypred/",eth[i],"/",trait[l],"/SBayesR.snpRes")
tar_con = 1
#load result
if(tar_coef%in%files==T){
  
  sbayes_result = fread(paste0(tar_coef),header = T)
  if(nrow(sbayes_result)>0){
     
    assoc = sbayes_result %>% 
      rename(BETA = A1Effect,
             SNP = Name) %>% 
      select(SNP, A1, BETA)
    write.table(assoc,file = paste0(temp.dir,"prs_prep"),col.names = T,row.names = F,quote=F)
    #calculate the PRS for SBayesR 
    res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                        "--score-col-nums 3 --threads 2 ",
                        "--score ",temp.dir,"prs_prep cols=+scoresums,-scoreavgs header no-mean-imputation ",
                        "--bfile ",ref_gene_pred,
                        " --out ",temp.dir,"PRS_SBayesR_tar"))
    
  }else{
    tar_con = 0
  }
}





#calculate the PRS for Polyfun
#if(i ==1){

poly_fun_file_out = paste0("/data/zhangh24/multi_ethnic/result/AOU/prs/polypred/",eth[1],"/",trait[l],"/poly_fun")
polyfun_result = fread(poly_fun_file_out)  

assoc = polyfun_result %>% 
  rename(BETA = BETA_MEAN) %>% 
  select(SNP, A1, BETA)
write.table(assoc,file = paste0(temp.dir,"prs_prep"),col.names = T,row.names = F,quote=F)
#calculate the PRS for SBayesR
res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                    "--score-col-nums 3 --threads 2 ",
                    "--score ",temp.dir,"prs_prep cols=+scoresums,-scoreavgs header no-mean-imputation ",
                    "--bfile ",ref_gene_pred,
                    " --out ",temp.dir,"PRS_polyfun"))

SBayesR_EUR = as.data.frame(fread(paste0(temp.dir, "PRS_SBayesR_EUR.sscore")))

polyfun_EUR = as.data.frame(fread(paste0(temp.dir, "PRS_polyfun.sscore")))
if(tar_con==1){
  SBayesR_tar = as.data.frame(fread(paste0(temp.dir, "PRS_SBayesR_tar.sscore")))
  prs_mat = cbind(SBayesR_EUR,polyfun_EUR[,5],SBayesR_tar[,5])
  colnames(prs_mat)[5:7] = c("SBayesR_EUR","Poly_fun","SBayesR_tar")
}else{
  
  prs_mat = cbind(SBayesR_EUR,polyfun_EUR[,5])
  colnames(prs_mat)[5:6] = c("SBayesR_EUR","Poly_fun")
}



prs_score = prs_mat[,5:ncol(prs_mat)]
colnames(prs_mat)[2] = "id"
pheno.dir = "/data/zhangh24/multi_ethnic/data/UKBB/phenotype/"
pheno_tuning = as.data.frame(fread(paste0(pheno.dir,trait[l],"/tuning+validation/",eth[i],"_tuning.txt")))
pheno_tuning = pheno_tuning[,1:2]
covar <- as.data.frame(fread(paste0(pheno.dir,"/covariates/tuning+validation/",eth[i],"_all_data.txt")))
pheno_tuning <- left_join(pheno_tuning, covar)
colnames(pheno_tuning) = c('id','y','sex','age',paste0('pc',1:10))
pheno_tuning_com = pheno_tuning[complete.cases(pheno_tuning$y),]
pheno_tuning = left_join(pheno_tuning_com,prs_mat,by = "id")
model.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_tuning)
y_tun = model.null$residual
prs_tun = pheno_tuning[,colnames(prs_score)]

pheno_vad = as.data.frame(fread(paste0(pheno.dir,trait[l],"/tuning+validation/",eth[i],"_validation.txt")))
pheno_vad = pheno_vad[,1:2]
pheno_vad <- left_join(pheno_vad, covar)
colnames(pheno_vad) = c('id','y','sex','age',paste0('pc',1:10))
pheno_vad_com = pheno_vad[complete.cases(pheno_vad$y),]
pheno_vad = left_join(pheno_vad_com,prs_mat,by = "id")
model.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_vad)
y_vad = model.null$residual
prs_vad = pheno_vad[,colnames(prs_score)]


  model = lm(y_tun~as.matrix(prs_tun))
  coef = coefficients(model)[-1] #take out intercept

prs = as.matrix(prs_vad)%*%coef

model = lm(y_vad~ prs)
r2 = summary(model)$r.square

data = data.frame(y = y_vad, x = prs)
R2Boot = function(data,indices){
  boot_data = data[indices, ]
  model = lm(y ~ x, data = boot_data)
  result = summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 = boot(data = data, statistic = R2Boot, R = 10000)

ci_result = boot.ci(boot_r2, type = "bca")


r2.result = data.frame(eth = eth[i],
                       trait = trait[l],
                       method = "PolyPred",
                       r2 = r2,
                       r2_low = ci_result$bca[4],
                       r2_high = ci_result$bca[5]
)
out.dir = paste0("/data/zhangh24/multi_ethnic/result/AOU/polypred/",eth[i],"/",trait[l],"/")
save(r2.result, file = paste0(out.dir, "polypred.result"))
#     system(paste0("rm -rf ", temp.dir))
#   }
# }
