#args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = 2
l = 1
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
system(paste0("cp /data/zhangh24/multi_ethnic/data/GBHS_plink/all_chr.bed ",temp.dir,"all_chr.bed"))
system(paste0("cp /data/zhangh24/multi_ethnic/data/GBHS_plink/all_chr.bim ",temp.dir,"all_chr.bim"))
system(paste0("cp /data/zhangh24/multi_ethnic/data/GBHS_plink/all_chr.fam ",temp.dir,"all_chr.fam"))

pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("overall","erpos","erneg")
setwd("/data/zhangh24/multi_ethnic/data/")
#load best TDLD result
method = "TDLD"
out.dir = paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/result/clump_result/")
load(paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/result/auc_tdld.rdata"))
result = auc.tdld[[2]]
idx <- which.max(result[,1])
r_ind = result$r2.vec[idx]
w_ind = result$wc.vec[idx]
k1 = result$p.eur.vec[idx]
k2 = result$p.tar.vec[idx]
method= "TDLD"
out.dir.prs = "/data/zhangh24/multi_ethnic/result/breast_cancer/prs/"
best.snp <- as.data.frame(fread(paste0(out.dir,"TDLD_rind_",r_ind,"_wcind_",w_ind,".clumped")))

r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)
out.dir = paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/result/clump_result/")
#load data from target ethnic group
load(paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS_mega.rdata"))
sum.tar = sum.data
sum.tar.select = sum.tar %>% 
  rename(SNP = V1,
         Effect_allele_TAR = Effect_allele,
         #Ref_allele_TAR = Alt_allele,
         BETA.TAR = Effect,
         SE.TAR = StdErr,
         BP = POS,
         ptar = P) %>% 
  select(CHR,SNP,BP,ID,BETA.TAR,SE.TAR,Effect_allele_TAR,
         ptar) 

#load data from EUR ethnic group
load("./AABC_data/BC_EUR_overall_mega_aligned.rdata")
sum.eur = sum.data
sum.eur.select = sum.eur %>% 
  rename(SNP=V1,
         peur = p_eur,
         BETA.EUR = Beta_eur_update,
         SE.EUR = Se_eur,
         Effect_allele_EUR = Effect_allele,
         #Ref_allele_EUR = Alt_allele,
         BP = POS) %>% 
  select(SNP,BETA.EUR,SE.EUR,Effect_allele_EUR,
         peur) 

sum.com <- left_join(sum.tar.select,sum.eur.select,by="SNP")
colSums(is.na(sum.com))

sum.com = sum.com %>% 
  mutate(BETA.EUR=ifelse(Effect_allele_TAR==Effect_allele_EUR,
                         BETA.EUR,
                         -BETA.EUR),
         Effect_allele_EUR =ifelse(Effect_allele_EUR==Effect_allele_TAR,
                        Effect_allele_EUR,Effect_allele_TAR)) %>% 
  mutate(z_stat_eur = BETA.EUR/SE.EUR,
         z_stat_tar = BETA.TAR/SE.TAR)
#align with best tdld prs
source("/data/zhangh24/multi_ethnic/code/stratch/EB_function.R")
best.prs.com = left_join(best.snp,sum.com)
beta_tar <- best.prs.com$BETA.TAR
sd_tar <- best.prs.com$SE.TAR
beta_eur <- best.prs.com$BETA.EUR
sd_eur <- best.prs.com$SE.EUR
EBprior = EstimatePrior(beta_tar,sd_tar,
                        beta_eur,sd_eur)
beta_tar <- sum.com$BETA.TAR
sd_tar <- sum.com$SE.TAR
beta_eur <- sum.com$BETA.EUR
sd_eur <- sum.com$SE.EUR
post_beta_mat = EBpost(beta_tar,sd_tar,beta_eur,sd_eur,EBprior)
post_beta_tar = post_beta_mat[,1,drop=F]
colnames(post_beta_tar) = "BETA"
sum.com$BETA = post_beta_tar


#get the min p-value for between the target ethnic group and EUR for shared snp
summary.com <-sum.com
temp = 1
out.dir.prs = "/data/zhangh24/multi_ethnic/result/breast_cancer/prs/"
n_pthres = length(pthres)
q_range = data.frame(
  paste0("p_value_",c(1:n_pthres)),
  rep(0,n_pthres),
  pthres,stringsAsFactors = F)
write.table(q_range,file = paste0(temp.dir,"2Dq_range_file"),row.names = F,col.names = F,quote=F)
out.dir.prs = "/data/zhangh24/multi_ethnic/result/breast_cancer/prs/"
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    LD <- as.data.frame(fread(paste0(out.dir,"TDLD_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    clump.snp <- LD
    prs.all <- left_join(clump.snp,summary.com)
    colSums(is.na(prs.all))
    prs.all = prs.all %>%
      rename(A1  = Effect_allele_TAR,
             P = ptar ) %>% 
      select(ID,CHR,BP,A1,BETA,P,peur) %>% 
      rename(SNP=ID)
    
    dup.id <- prs.all$SNP[duplicated(prs.all$SNP)]
    if(length(dup.id)!=0){
      #keep more space in case of multiple duplication
      remove.idx = rep(0,2*length(dup.id))
      temp = 1
      for(k in 1:length(dup.id)){
        jdx <- which(prs.all$SNP==dup.id[k])
        dup.length = length(jdx[-which.min(prs.all$P[jdx])])
        remove.idx[temp:(temp-1+dup.length)] = jdx[-which.min(prs.all$P[jdx])]
        temp = temp+dup.length
      }
      remove.idx = remove.idx[1:(temp-dup.length)]
    }
    prs.all = prs.all[-remove.idx,]
    
    
    for(k1 in 1:length(pthres)){
      
      #keep al the SNPs with peur pass the threshold
      prs.file <- prs.all %>% 
        mutate(P = replace(P,peur<=pthres[k1],1E-20))%>%
        select(SNP,A1,BETA,P)
      write.table(prs.file,file = paste0(temp.dir,"2DEBprs_coeff"),col.names = T,row.names = F,quote=F)
      
      p.value.file = prs.file %>% 
        select(SNP,P)
      
      
      write.table(p.value.file,file = paste0(temp.dir,"2DEBp_value"),col.names = T,row.names = F,quote=F)
      # # com.prs = left_join(prs.file,p.value.file,by="SNP")
      # com.prs.filter = com.prs %>%
      #   filter(P<=pthres[9])
      #bim <- fread(paste0(data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",j,".bim"))
      # idx <- which(bim$V2=="rs1967017")
      # bim[idx,]
      if(nrow(prs.file)>0){
        res <- system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",temp.dir,"all_chr --q-score-range ",temp.dir,"2Dq_range_file ",temp.dir,"2DEBp_value header --score ",temp.dir,"2DEBprs_coeff header no-sum no-mean-imputation  --out ",temp.dir,"prs_EB",trait[l],"_rind_",r_ind,"_wcind_",w_ind,"p_value_",k1))
        #res <- system(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic_data_analysis/plink2 --q-score-range ",temp.dir,"q_range_file ",temp.dir,"p_value_chr_",j," header --threads 2 --score ",temp.dir,"prs_coeff_chr_",j," header no-mean-imputation --bfile ",data.dir,trait[1],"/",eth[i],"/geno/mega/chr.qc",j," --out ",temp.dir,"prs_chr_",j))
        print("step2 finished")
        res = system(paste0("mv ",temp.dir,"*.profile ",out.dir.prs))
        #system(paste0("/data/zhangh24/software/plink2 --score ",cur.dir,eth[i],"/prs/prs_file_pvalue_",k,"_rho_",l,"_size_",m,,"_rep_",i_rep," no-sum no-mean-imputation --bfile ",cur.dir,eth[i],"/all_chr.tag --exclude /data/zhangh24/multi_ethnic/result/LD_simulation/",eth[i],"/duplicated.id  --out ",cur.dir,eth[i],"/prs/prs_",k,"_rho_",l,"_size_",m))
        if(res==2){
          stop()
        }
        
      }
    }
  }
}