args = commandArgs(trailingOnly = T)
#l =1,2,3 for overall, ERpos, ERneg using BCAC
#l = 4 for overall using BCAC and UKB meta
l = as.numeric(args[[1]])
#best eur prs
AUCEstimate <- function(data, indices,tar_trait,scorename) {
  d <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status=tar_trait, variable=scorename,
                        confounders=~EV1+EV2+EV3+EV4+EV5+
                          EV6+EV7+EV8+EV9+EV10,
                        data=d, precision=seq(0.1,0.9, by=0.1))
  return(roc_obj$auc)
} 

setwd("/data/zhangh24/multi_ethnic/data/")
library(data.table)
library(tidyverse)
library(RISCA)
library(boot)

#load GHBS_id
bim <- fread("/data/zhangh24/multi_ethnic/data/GBHS_plink/all_chr.bim",header=F) 
bim = bim %>% 
  unite("chr.pos",V1,V4,sep=":",remove=F)
bim = bim %>% 
  rename(GBHS_id = V2)


trait = c("overall","erpos","erneg","overall_ukb_meta")

prs.snp <- read.csv("/data/zhangh24/multi_ethnic/data/breast_cancer_330SNPs.csv",header=T)
prs.snp = prs.snp %>% 
  unite("chr.pos",CHR,Position,sep=":",remove=F)

prs.snp = prs.snp %>% 
  select(variant,chr.pos,Reference.allele,Effect.allele)

load(paste0("./AABC_data/BC_EUR_",trait[l],"_aligned.rdata"))
#match best eur SNP with AABC summary level statistics
prs.snp.eur = left_join(prs.snp,sum.data,by="chr.pos") %>% 
  filter(((Effect_allele==Effect.allele)&(Alt_allele==Reference.allele))|
           (Effect_allele==Reference.allele)&(Alt_allele==Effect.allele)) %>% 
  select(c(colnames(sum.data),"chr.pos","variant"))

# jdx <- which(bim$V4==106358620)
# idx <- which(sum.data$POS==187503758 )
#idx <- which(prs.snp.eur.temp$ID%in%prs.snp.eur$ID==F)
#idx <- which(sum.data$chr.pos =="4:84370124")
#idx <- which(sum.data.temp$chr.pos =="4:84370124")
#match best eur SNP with AABC summary level statistics and Ghana genotype data
prs.snp.eur.update = left_join(prs.snp.eur,bim,by="chr.pos") %>% 
  filter(((Effect_allele==V5)&(Alt_allele==V6))|
           (Effect_allele==V6)&(Alt_allele==V5))
#prepare eur coefficients
prs.file = prs.snp.eur.update %>% 
  select(GBHS_id,Effect_allele,Beta_eur_update) %>% 
  rename(SNP=GBHS_id,
         A1 = Effect_allele,
         BETA = Beta_eur_update)
prs.eur = prs.file

#prepare afr coefficients
if(l==4){
  load(paste0("./AABC_data/BC_AFR_",trait[1],"remove_GHBS.rdata"))  
}else{
  load(paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS.rdata"))  
}


# sum.data = sum.data %>% 
#   unite("chr.pos",CHR,POS,sep=":",remove=F)


prs.snp.tar = left_join(prs.snp,sum.data,by="chr.pos") %>% 
  filter(((Effect_allele==Effect.allele)&(Alt_allele==Reference.allele))|
           (Effect_allele==Reference.allele)&(Alt_allele==Effect.allele)) %>% 
  select(c(colnames(sum.data),"chr.pos","variant"))


prs.snp.tar.update = left_join(prs.snp.tar,bim,by="chr.pos") %>% 
  filter(((Effect_allele==V5)&(Alt_allele==V6))|
           (Effect_allele==V6)&(Alt_allele==V5))
prs.file = prs.snp.tar.update %>% 
  select(GBHS_id,Effect_allele,Effect) %>% 
  rename(SNP=GBHS_id,
         A1 = Effect_allele,
         BETA = Effect)
prs.tar = prs.file

all.equal(prs.tar$SNP,prs.eur$SNP)
all.equal(prs.tar$A1,prs.eur$A1)


#prepare EB coefficients
source("/data/zhangh24/multi_ethnic/code/stratch/EB_function.R")

beta_tar <- prs.tar$BETA
sd_tar <- prs.snp.tar.update$StdErr
beta_eur <- prs.eur$BETA
sd_eur <- prs.snp.eur.update$Se_eur
EBprior = EstimatePrior(beta_tar,sd_tar,
                        beta_eur,sd_eur)

post_beta_mat = EBpost(beta_tar,sd_tar,beta_eur,sd_eur,EBprior)
post_beta_tar = post_beta_mat[,1,drop=F]

BETA_EB = post_beta_tar
#prs.tar.all$BETA = post_beta_tar

prs.file =cbind(prs.eur,BETA_tar=prs.tar$BETA,BETA_EB)

# idx <- which(prs.snp.tar$chr.pos%in%prs.snp.eur$chr.pos==F)
# jdx = which(sum.data$chr.pos=="7:74341926")

write.table(prs.file,file = paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/prs/best_eurprs_",trait[l]),row.names = F,col.names = T,quote=F)
out.dir.prs = "/data/zhangh24/multi_ethnic/result/breast_cancer/prs/"
data.dir = "/data/zhangh24/multi_ethnic/data/"
res <- system(paste0("/data/zhangh24/software/plink2_alpha --score-col-nums 3,4,5 --score ",out.dir.prs,"best_eurprs_",trait[l]," header no-mean-imputation --bfile ",data.dir,"GBHS_plink/all_chr --out ",out.dir.prs,"prs_best_eur_",trait[l]))

#calculate AUC for different traits
if(l==1){
auc.est = rep(0,3)
auc.se = rep(0,3)
names(auc.est) = c("EUR_coef","tar_coef","EB_coef")
for(k in 1:3){
  prs.score = fread(paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/prs/prs_best_eur_",trait[l],".sscore"))
  prs.score = prs.score %>% 
    separate(IID,into=c("ID","ohter_id","nci_id"),remove=F)
  
  prs.score = prs.score %>% select(ID,paste0("SCORE",k,"_AVG"))
  pheno <- fread(paste0(data.dir,"GBHS_pheno.csv"))
  colnames(pheno)[5] = "ER_status"
  
    pheno = pheno %>% 
      mutate(overall_status=ifelse(Status=="Control",0,1))
    
    pheno.update = left_join(pheno,prs.score,by=c("SB_ID"="ID"))
    roc_obj <- roc.binary(status="overall_status", variable=paste0("SCORE",k,"_AVG"),
                          confounders=~EV1+EV2+EV3+EV4+EV5+
                            EV6+EV7+EV8+EV9+EV10,
                          data=pheno.update, precision=seq(0.1,0.9, by=0.1))
    boot_result <- boot(data = pheno.update,statistic = AUCEstimate,
                        R = 3000,scorename = paste0("SCORE",k,"_AVG"),
                        tar_trait = "overall_status")
    boot.ci(boot_result, type="bca")
    #roc_obj = roc(pheno.update$ERneg,pheno.update$SCORE)
    auc.est[k] = roc_obj$auc
    auc.se[k] = sqrt(var(boot_result[[2]]))
    
  }
}else if(l==2){
    auc.est = rep(0,3)
    auc.se = rep(0,3)
    names(auc.est) = c("EUR_coef","tar_coef","EB_coef")
    for(k in 1:3){
      prs.score = fread(paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/prs/prs_best_eur_",trait[l],".sscore"))
      prs.score = prs.score %>% 
        separate(IID,into=c("ID","ohter_id","nci_id"),remove=F)
      
      prs.score = prs.score %>% select(ID,paste0("SCORE",k,"_AVG"))
      pheno <- fread(paste0(data.dir,"GBHS_pheno.csv"))
      colnames(pheno)[5] = "ER_status"
      
    pheno = pheno %>% 
      mutate(ERpos=case_when(Status == "Control" ~ 0,
                             ER_status== 1 ~ 1,
                             TRUE ~ 888)) %>% 
      filter(ERpos!=888)
    
    pheno.update = left_join(pheno,prs.score,by=c("SB_ID"="ID"))
    roc_obj <- roc.binary(status="ERpos", variable=paste0("SCORE",k,"_AVG"),
                          confounders=~EV1+EV2+EV3+EV4+EV5+
                            EV6+EV7+EV8+EV9+EV10,
                          data=pheno.update, precision=seq(0.1,0.9, by=0.1))
    boot_result <- boot(data = pheno.update,statistic = AUCEstimate,
                        R = 3000,scorename = paste0("SCORE",k,"_AVG"),
                        tar_trait = "ERpos")
    boot.ci(boot_result, type="bca")
    auc.est[k] = roc_obj$auc
    auc.se[k] = sqrt(var(boot_result[[2]]))
  }
  
  }else if(l==3){
    auc.est = rep(0,3)
    auc.se = rep(0,3)
    names(auc.est) = c("EUR_coef","tar_coef","EB_coef")
    for(k in 1:3){
      prs.score = fread(paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/prs/prs_best_eur_",trait[l],".sscore"))
      prs.score = prs.score %>% 
        separate(IID,into=c("ID","ohter_id","nci_id"),remove=F)
      
      prs.score = prs.score %>% select(ID,paste0("SCORE",k,"_AVG"))
      pheno <- fread(paste0(data.dir,"GBHS_pheno.csv"))
      colnames(pheno)[5] = "ER_status"
      
      pheno = pheno %>% 
        mutate(ERneg=case_when(Status == "Control" ~ 0,
                               ER_status== 2 ~ 1,
                               TRUE ~ 888)) %>% 
        filter(ERneg!=888)
      
      
      
      pheno.update = left_join(pheno,prs.score,by=c("SB_ID"="ID"))
      
     
      
      roc_obj <- roc.binary(status="ERneg", variable=paste0("SCORE",k,"_AVG"),
                         confounders=~EV1+EV2+EV3+EV4+EV5+
                           EV6+EV7+EV8+EV9+EV10,
                         data=pheno.update, precision=seq(0.1,0.9, by=0.1))
      boot_result <- boot(data = pheno.update,statistic = AUCEstimate,
                          R = 3000,scorename = paste0("SCORE",k,"_AVG"),
                          tar_trait = "ERneg")
      boot.ci(boot_result, type="bca")
      #roc_obj = roc(pheno.update$ERneg,pheno.update$SCORE)
      auc.est[k] = roc_obj$auc
      auc.se[k] = sqrt(var(boot_result[[2]]))
      
    }
  }else if(l==4){
    
    auc.est = rep(0,3)
    auc.se = rep(0,3)
    names(auc.est) = c("EUR_coef","tar_coef","EB_coef")
    for(k in 1:3){
      prs.score = fread(paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/prs/prs_best_eur_",trait[l],".sscore"))
      prs.score = prs.score %>% 
        separate(IID,into=c("ID","ohter_id","nci_id"),remove=F)
      
      prs.score = prs.score %>% select(ID,paste0("SCORE",k,"_AVG"))
      pheno <- fread(paste0(data.dir,"GBHS_pheno.csv"))
      colnames(pheno)[5] = "ER_status"
      
      pheno = pheno %>% 
        mutate(overall_status=ifelse(Status=="Control",0,1))
      
      pheno.update = left_join(pheno,prs.score,by=c("SB_ID"="ID"))
      roc_obj <- roc.binary(status="overall_status", variable=paste0("SCORE",k,"_AVG"),
                            confounders=~EV1+EV2+EV3+EV4+EV5+
                              EV6+EV7+EV8+EV9+EV10,
                            data=pheno.update, precision=seq(0.1,0.9, by=0.1))
      boot_result <- boot(data = pheno.update,statistic = AUCEstimate,
                          R = 3000,scorename = paste0("SCORE",k,"_AVG"),
                          tar_trait = "overall_status")
      boot.ci(boot_result, type="bca")
      #roc_obj = roc(pheno.update$ERneg,pheno.update$SCORE)
      auc.est[k] = roc_obj$auc
      auc.se[k] = sqrt(var(boot_result[[2]]))
      
    }
  }

auc.result = list(auc.est,auc.se)
save(auc.result,file = paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/result/auc_result_",l))












# attributes(temp_obj)
# n_p = sum(pheno.update$overall_status)
# n_n = sum(1-pheno.update$overall_status)
# library(auctestr)
# se_auc(as.numeric(auc(roc_obj)),n_p,n_n)
# temp_obj$se
# ci(roc_obj)
# ci.se(roc_obj)
# idx <- which(prs.snp.eur$ID%in%prs.snp.eur.update$ID==F)
# prs.snp.eur[idx,]
