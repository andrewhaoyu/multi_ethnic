#goal: preprare the prs for 23andMe analyses
#prs-csx five ancestries
#xpass
#polypred


args = commandArgs(trailingOnly = T)
#i represent ethnic group
#l represent chromosome

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)

eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")
#load snp_info file to get corresponding im.id
data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"
load(paste0(data.dir,"snpinfo/snpinfo_mega.RData"))
snpinfo_mega_filter = snpinfo_mega %>% 
  filter(!is.na(im.data.id)) %>% 
  select(im.data.id,assay.name)
#######load all the methods and find the unique SNP ID
sum = as.data.frame(fread(paste0(data.dir,eth[i],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
#use orignal summary stat A1 as our effect allele
sum_infor = sum %>% 
  select(rsid, A1, A2) %>% 
  rename(SNP = rsid)
#load xpass results

method = "XPASS"
out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/",method,"/",eth[i],"/",trait[l],"/")
file_out = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/XPASS/",eth[i],"/",trait[l],"/")
load(paste0(file_out, "_param.RData"))
SNP_XPASS = fit_bbj$mu[,"SNP"]

#load Polypred results: SBAYESR EUR, SBAYESR target, Polyfun EUR

method = "Polypred"
result = fread(paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/polypred/",eth[1],"/",trait[l],"/SBayesR.snpRes"))
SNP_SBAYES_EUR = result$Name
result = fread(paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/polypred/",eth[i],"/",trait[l],"/SBayesR.snpRes"))
SNP_SBAYES_tar = result$Name
result = fread(paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/polypred/",eth[1],"/",trait[l],"/poly_fun"))
SNP_polyfun = result$SNP

#load prs-csx five ancestries results
out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/PRSCSx/EUR/",trait[l],"/")
#   for(l in 1:7){
phi = c("1e+00","1e-02","1e-04","1e-06")
v = 1
snp.id.list = list()
for(k in 1:5){
  load(paste0(out.dir.prs,"sum_five_",eth[k],"_pst_eff_a1_b0.5_phi",phi[v],".rdata"))
  snp.id.list[[i]] = data.frame(SNP = data$V2)
  
}

SNP_PRS_csx = rbindlist(snp.id.list)


#combine all SNPs together
SNP = data.frame(SNP = unique(c(SNP_XPASS,
               SNP_SBAYES_EUR,
               SNP_SBAYES_tar,
               SNP_polyfun,
               SNP_PRS_csx$SNP)))
#only keep SNPs that exist in the target population
SNP_infor = inner_join(SNP, sum_infor, by = "SNP")
snp_id = left_join(SNP_infor,
                   snpinfo_mega_filter,
                   by=c("SNP"="assay.name"))
#flip the allele to be the right order as the original summary statistics

#load xpass results

method = "XPASS"
out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/",method,"/",eth[i],"/",trait[l],"/")
file_out = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/XPASS/",eth[i],"/",trait[l],"/")
load(paste0(file_out, "_param.RData"))
BETA = fit_bbj$mu[,c("SNP","A1","mu_XPASS1")]
colnames(BETA) = c("SNP", "effect_allele", "beta")
# idx <- which(temp_data$effect_allele!=A1)
# jdx = which(temp_data$effect_allele[idx]!=temp_data$A2[idx])
# temp_data[idx[jdx[1:10]],]
AlignSNP = function(snp_id,BETA){
  temp_data = left_join(snp_id, BETA, by = c("SNP" = "SNP"))
  temp_data = temp_data %>% 
    mutate(beta_update = case_when(
      is.na(effect_allele) == T ~ 0,
      effect_allele == A2 ~ -beta,
      effect_allele == A1 ~ beta,
      effect_allele != A1 & effect_allele != A2 ~ 0 #insertation and deletion, not match with reference 
    ))
return(temp_data)  
}

temp_data = AlignSNP(snp_id, BETA)
beta_xpass = temp_data$beta_update


#load Polypred results: SBAYESR EUR, SBAYESR target, Polyfun EUR

method = "Polypred"
result = fread(paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/polypred/",eth[1],"/",trait[l],"/SBayesR.snpRes"))
BETA = result[,c("Name", "A1", "A1Effect")]
colnames(BETA) = c("SNP", "effect_allele", "beta")
temp_data = AlignSNP(snp_id,BETA)
beta_sbayes_eur = temp_data$beta_update
result = fread(paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/polypred/",eth[i],"/",trait[l],"/SBayesR.snpRes"))
BETA = result[,c("Name", "A1", "A1Effect")]
colnames(BETA) = c("SNP", "effect_allele", "beta")
temp_data = AlignSNP(snp_id,BETA)
beta_sbayes_tar = temp_data$beta_update
result = fread(paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/polypred/",eth[1],"/",trait[l],"/poly_fun"))
SNP_polyfun = result[,c("SNP", "A1", "BETA_MEAN")]
colnames(BETA) = c("SNP", "effect_allele", "beta")
temp_data = AlignSNP(snp_id,BETA)
beta_sbayes_polyfun = temp_data$beta_update

beta_polypred = cbind(beta_sbayes_eur,
                      beta_sbayes_tar,
                      beta_sbayes_polyfun)



out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/PRSCSx/EUR/",trait[l],"/")
BETA_mat = matrix(0, nrow = nrow(snp_id),ncol = length(eth)*length(phi))
#   for(l in 1:7){
temp = 1
for(v in 1:length(phi))
for(k in 1:length(eth)){
  
  load(paste0(out.dir.prs,"sum_five_",eth[k],"_pst_eff_a1_b0.5_phi",phi[v],".rdata"))
  BETA = data[,c("V2", "V4", "V6")]
  colnames(BETA) = c("SNP", "effect_allele", "beta")
  temp_data = AlignSNP(snp_id,BETA)
  BETA_mat[,temp] = temp_data$beta_update
  temp = temp + 1
  
  
  
}



data = cbind(snp_id,beta_xpass,beta_polypred,BETA_mat)

prs.snp = data[,c(4,5:ncol(data))]
out.dir.organize.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/organize_prs/three_method/",eth[i],"/",trait[l],"/")
write.table(prs.snp,file = paste0(out.dir.organize.prs,"prs.file"),row.names = F,col.names = F,quote=F)
#   }


# 
# temp = fread( paste0(out.dir.organize.prs,"prs.file"))
# # }
# 







for(i in 1:5){
  for(l in 1:7){
    file_dir = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/polypred/",eth[i],"/",trait[l],"/")
    files = dir(file_dir, pattern = "SBayesR.snpRes")
    if(length(files)==0){
      print(c(i,l))
    }
    
  }
}
result = fread(paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/polypred/",eth[i],"/",trait[l],"/SBayesR.snpRes"))