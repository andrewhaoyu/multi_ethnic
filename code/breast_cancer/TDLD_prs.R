args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = 2
l = 1
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("overall","erpos","erneg")
setwd("/data/zhangh24/multi_ethnic/")
r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)
method = "TDLD"
out.dir = paste0("/data/zhangh24/multi_ethnic/result/breast_cancer/result/clump_result/")
load("./AABC_data/BC_EUR_overall_mega_aligned.rdata")
sum.eur = sum.data
sum.eur.select = sum.eur %>% 
  rename(SNP=V1,
         peur = p_eur) %>% 
  select(SNP,peur) 
load(paste0("./AABC_data/BC_AFR_",trait[l],"remove_GHBS_mega.rdata"))
sum.tar = sum.data
sum.data.assoc = sum.tar %>% 
  rename(SNP = V1,
         A1 = Effect_allele,
         BETA = Effect,
         BP = POS) %>% 
  select(CHR,SNP,BP,A1,BETA,P) 
#get the min p-value for between the target ethnic group and EUR for shared snp
summary.com <- left_join(sum.data.assoc,sum.eur.select,by="SNP")
temp = 1
out.dir.prs = "/data/zhangh24/multi_ethnic/result/breast_cancer/prs/"
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    LD <- as.data.frame(fread(paste0(out.dir,"TDLD_rind_",r_ind,"_wcind_",w_ind,".clumped")))
    clump.snp <- LD
    prs.all <- left_join(clump.snp,summary.com)
    
    
    for(k1 in 1:length(pthres)){
      for(k2 in 1:length(pthres)){
        prs.temp <- prs.all %>% 
          filter((P<=pthres[k1]|
                    peur<=pthres[k2])) %>% 
          select(SNP,BETA)
        if(nrow(prs.temp)>0){
          write.table(prs.temp,file = paste0(out.dir.prs,temp,"_",method,"_rind_",r_ind,"_wcind_",w_ind,"_ptar_",k1,"_peur_",k2),row.names = F,col.names = T,quote=F)  
          temp = temp + 1
        }
        
      }
      
    }
  }
}