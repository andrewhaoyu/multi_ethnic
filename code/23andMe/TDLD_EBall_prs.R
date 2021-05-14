args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
#generate index table
method = "tdld"
#load all results
library(data.table)
library(tidyverse)
out.dir = "/data/zhangh24/multi_ethnic/result/cleaned/result_summary/"
data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"
eth_group = c("european","african_american",
              "latino","east_asian","south_asian")

eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")

r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
temp = 1
total = length(r2_vec)*length(wc_base_vec)*length(pthres)*length(pthres)
r2 = rep(0,total)
wc = rep(0,total)
ptar = rep(0,total)
peur = rep(0,total)



for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    for(k1 in 1:length(pthres)){
      for(k2 in 1:length(pthres)){
        r2[temp] = r_ind
        wc[temp] = w_ind
        ptar[temp] = k1
        peur[temp] = k2
        temp = temp+1
      }
    }
  }
}

index.table = data.frame(r2,wc,ptar,peur)
#load best TDLD result
method = "tdld"
result <- read.csv(paste0(out.dir,eth_group[i],"_",trait[l],"_",method))
idx <- index.table[which.max(result),]
temp = which.max(result)
r_ind = idx[1]
w_ind = idx[2]
k1 = idx[3]
k2 = idx[4]
method= "TDLD"
out.dir.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/",method,"/",eth[i],"/",trait[l],"/")
best.prs = read.table(paste0(out.dir.prs,temp,"_",method,"_rind_",r_ind,"_wcind_",w_ind,"_ptar_",k1,"_peur_",k2),header=T)
#load data from target ethnic group
sum.tar = as.data.frame(fread(paste0(data.dir,eth[i],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
sum.tar.select = sum.tar %>% 
  select(rsid,BETA,SD,P) %>% 
  rename(BETA.TAR = BETA,
         SD.TAR = SD)
#load data from EUR ethnic group
sum.eur = as.data.frame(fread(paste0(data.dir,eth[1],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
sum.eur.select = sum.eur %>% 
  select(rsid,BETA,SD,P) %>% 
  rename(BETA.EUR= BETA,
         SD.EUR = SD,
         peur = P)
#combine the two ethnic group data
summary.com = left_join(sum.tar.select,sum.eur.select,by = "rsid")
#load other ethnic group data
i_can <- setdiff(c(2:5),i)
beta.mat.list = list()
sd.mat.list = list()
temp = 1
for(i_sub in i_can){
  sum.data <- as.data.frame(fread(paste0(data.dir,eth[i_sub],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
  colnames(sum.data)[2] <- "SNP"
  sum.data.select = sum.data %>% 
    select(rsid,BETA,SD) %>% 
    rename(BETA.sub = BETA,
           SD.sub = SD)
  summary.com.temp <- left_join(summary.com,sum.data.select,by="rsid")
  # #match allele
  # summary.com.temp = summary.com.temp %>% 
  #   mutate(BETA.sub.new = ifelse(A1==A1.sub,BETA.sub,-BETA.sub)) 
  beta.mat.list[[temp]] =  summary.com.temp %>% select(BETA.sub)
  sd.mat.list[[temp]] = summary.com.temp %>% select(SD.sub)
  temp = temp+1
}
beta.mat = bind_cols(beta.mat.list)
sd.mat = bind_cols(sd.mat.list)

summary.com$beta.mat = beta.mat
summary.com$sd.mat = sd.mat


#
best.prs.com = left_join(best.prs,summary.com,by=c("SNP"="rsid"))

beta_tar <- best.prs.com$BETA.TAR
sd_tar <- best.prs.com$SD.TAR
beta_eur <- best.prs.com$BETA.EUR
sd_eur <- best.prs.com$SD.EUR
beta_other_mat = best.prs.com$beta.mat
sd_other_mat = best.prs.com$sd.mat
source("/data/zhangh24/multi_ethnic/code/stratch/EB_function.R")
EBPrior = EstimatePriorMulti(beta_tar,sd_tar,
                                 beta_eur,sd_eur,
                                 beta_other_mat,
                                 sd_other_mat)
beta_tar <- summary.com$BETA.TAR
sd_tar <- summary.com$SD.TAR
beta_eur <- summary.com$BETA.EUR
sd_eur <- summary.com$SD.EUR
beta_other_mat = summary.com$beta.mat
sd_other_mat = summary.com$sd.mat


post_beta_mat = EBpostMulti(beta_tar,sd_tar,
                            beta_eur,sd_eur,beta_other_mat,
                            sd_other_mat,
                            EBPrior)
post_beta_tar = post_beta_mat[,1,drop=F]



colnames(post_beta_tar) = "BETA"
summary.com$BETA = post_beta_tar
colnames(summary.com)[1] = "SNP"

#install.packages("dplyr")
#install.packages("vctrs")
#for(i in 1:length(eth)){
# for(l in 1:length(trait)){
out.dir = paste0("/data/zhangh24/multi_ethnic/result/cleaned/clumping_result/TDLD/",eth[i],"/",trait[l],"/")
temp = 1
method = "TDLD_EBall"
out.dir.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/",method,"/",eth[i],"/",trait[l],"/")
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
# }
#}
