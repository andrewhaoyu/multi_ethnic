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

# for(i in 1:2){
#   for(l in 1:3){
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/GLGC_cleaned/"
out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/",eth,"/",trait,"/")
#########LD clumping#################
#load gwas summary statistics
sum.data = as.data.frame(fread(paste0(data.dir,eth,"/",trait,".txt"),header=T))
# write.table(sum.data.MAF,file = paste0("/lscratch/",sid,"/test/",eth[i],"_summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".out")
#             ,col.names = T,row.names = F,quote=F) 
#prepare association file for plink
sum.data.assoc = sum.data %>% 
  mutate(P = as.numeric(P)) %>% 
  rename(SNP = rsID,
         BP = POS_b37) %>% 
  select(CHR,SNP,BP,A1,BETA,P) 

#sum.data.assoc = sum.data.assoc[,c("CHR","SNP","BP","A1","BETA","P")]
#idx <- which(sum.data.assoc$SNP=="rs4970836")
write.table(sum.data.assoc,file = paste0(temp.dir,"all_chr_assoc.txt"),col.names = T,row.names = F,quote=F)

# dim(summary)
# head(summary)
pthr = 1
r2thr = 0.1
kbpthr = 500#cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
#code <- rep("c",5*3*3)
#system(paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/",eth[i],"/chr_all --clump ",cur.dir,eth[i],"/summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,".out --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",cur.dir,eth[i],"/LD_clump_rho_",l,"_size_",m,"_rep_",i_rep))
#plink2 is actually PLINK1.9
res = system(paste0("/data/zhangh24/software/plink2 --bfile ",temp.dir,eth,"all_chr --clump ",temp.dir,"all_chr_assoc.txt --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",temp.dir,"LD_clump"))
system(paste0("mv ",temp.dir,"LD_clump.clumped ",out.dir))
if(res==2){
  stop()
}
print("LD clumping end")
#########LD clumping ended#################

#########PRS calculation#################
temp.dir = paste0('/lscratch/',sid,'/test/')
system(paste0("mkdir ",temp.dir,"ukb"))
geno.data = paste0("/data/zhangh24/multi_ethnic/data/UKBB/genotype/all_data/")
system(paste0("cp ",geno.data,eth,"/all_chr.bed ",temp.dir,"ukb/all_chr.bed"))
system(paste0("cp ",geno.data,eth,"/all_chr.bim ",temp.dir,"ukb/all_chr.bim"))
system(paste0("cp ",geno.data,eth,"/all_chr.fam ",temp.dir,"ukb/all_chr.fam"))


pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
LD <- as.data.frame(fread(paste0(out.dir,"LD_clump.clumped")))
clump.snp <- LD[,3,drop=F]
prs.all <- left_join(clump.snp,sum.data.assoc)
temp = 1
out.dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/clumping_result/PT/",eth,"/",trait,"/")
out.dir.prs <- paste0("/data/zhangh24/multi_ethnic/result/GLGC/prs/PT/",eth,"/",trait,"/")


n_pthres <- length(pthres)
#q_range = data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres),stringsAsFactors = F)

prs.file = prs.all[,c("SNP","A1","BETA")]

write.table(prs.file,file = paste0(temp.dir,"prs_coeff"),col.names = T,row.names = F,quote=F)
p.value.file = prs.all[,c("SNP","P")]

write.table(p.value.file,file = paste0(temp.dir,"p_value"),col.names = T,row.names = F,quote=F)

q_range = data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres))

temp = 1
for(k in 1:length(pthres)){
  q_range[temp,1] = paste0("p_value_",k)
    q_range[temp,3] = pthres[k]
    temp = temp + 1
}
q_range = q_range[1:(temp-1),]
write.table(q_range,file = paste0(temp.dir,"q_range_file"),row.names = F,col.names = F,quote=F)

#plink2_alpha is plink2
#PRS = G*beta/(2*number of SNPs) #column header is SCORE_AVG
#PRS = G*beta
res <- system(paste0("/data/zhangh24/software/plink2_alpha --q-score-range ",temp.dir,"q_range_file ",temp.dir,"p_value header --threads 2 --score ",temp.dir,"prs_coeff cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile ",temp.dir,"ukb/all_chr --out ",temp.dir,"prs"))
system(paste0("ls ",temp.dir,""))
if(res==2){
    stop()
}
print("PRS calculation ended")



prs_list = list()
temp = 1
for(k in 1:length(pthres)){
  
    prs_temp = fread(paste0(temp.dir,"prs.p_value_",k,".sscore"))
    # times (2*number of SNPs)
    prs_list[[temp]] = prs_temp[,5:ncol(prs_temp)]
    
    colnames(prs_list[[temp]]) = paste0("p_value_",k)
    temp = temp + 1
  
}
prs_mat = as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
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
r2_tun_vec = rep(0,length(pthres))
model.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_tuning)
for(k in 1:length(pthres)){
  prs = pheno_tuning[,paste0("p_value_",k)]
  model.prs <- lm(model.null$residual~prs,data=pheno_tuning)
  r2_tun_vec[k] = summary(model.prs)$r.square
}
#find best idx
idx = which.max(r2_tun_vec)
#evaluate on validation
model.vad.null  =  lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_vad)
prs = pheno_vad[,paste0("p_value_",idx)]
model.vad.prs <- lm(model.vad.null$residual~prs,data=pheno_vad)
r2 = summary(model.vad.prs)$r.square


r2.result = data.frame(eth = eth,
                       trait = trait,
                       method = "CT",
                       r2 = r2
)


ct.result = list(r2.result,r2_tun_vec)


save(ct.result, file = paste0(out.dir, "CT.result"))
#find best cutoff for EUR by using all data as tuning


#########R2 calculation################# 
# 
# 
pheno.dir = "/data/zhangh24/multi_ethnic/data/UKBB/phenotype/"
pheno = as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",eth,"_all_data.txt")))
pheno = pheno[,1:2]
covar <- as.data.frame(fread(paste0(pheno.dir,"/covariates/tuning+validation/",eth,"_all_data.txt")))
pheno <- left_join(pheno, covar)
colnames(pheno) = c('id','y','sex','age',paste0('pc',1:10))
pheno_com = pheno[complete.cases(pheno$y),]

#combine prs with pheno_all
pheno_all = left_join(pheno_com,prs_mat,by = "id")

n_fold = 20
r2_vec = rep(0,n_fold)
for (fold in 1:n_fold){
  set.seed(123*fold)
  ids1 = sample(1:nrow(pheno_all), ceiling(nrow(pheno_all)/2), replace = F)
  ids2 = setdiff(1:nrow(pheno_all), ids1)

  pheno_tuning = pheno_all[ids1,]
  pheno_validation = pheno_all[ids2,]
  r2_tun_vec = rep(0,length(pthres))
  #calculate R2 for each of the tuning dataset
  model.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_tuning)
  for(k in 1:length(pthres)){
    prs = pheno_tuning[,paste0("p_value_",k)]
    model.prs <- lm(model.null$residual~prs,data=pheno_tuning)
    r2_tun_vec[k] = summary(model.prs)$r.square
  }
  #find best idx
  idx = which.max(r2_tun_vec)
  #evaluate on validation
  model.vad.null  =  lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_validation)
  prs = pheno_validation[,paste0("p_value_",idx)]
  model.vad.prs <- lm(model.vad.null$residual~prs,data=pheno_validation)
  r2_vec[fold] = summary(model.vad.prs)$r.square

}
r2 = mean(r2_vec)

r2.result = data.frame(eth = eth,
                       trait = trait,
                       method = "CT",
                       r2 = r2
)
#find best cutoff for EUR by using all data as tuning

  pheno_tuning = pheno_all
  r2_tun_vec = rep(0,length(pthres))
  #calculate R2 for each of the tuning dataset
  model.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_tuning)
  for(k in 1:length(pthres)){
    prs = pheno_tuning[,paste0("p_value_",k)]
    model.prs <- lm(model.null$residual~prs,data=pheno_tuning)
    r2_tun_vec[k] = summary(model.prs)$r.square
  }

  ct.result = list(r2.result,r2_tun_vec)


save(ct.result, file = paste0(out.dir, "CT.result"))
#   }
# }