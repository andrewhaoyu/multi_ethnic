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

#load xpass results

method = "XPASS"
out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/",method,"/",eth[i],"/",trait[l],"/")
file_out = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/XPASS/",eth[i],"/",trait[l],"/")
load(paste0(file_out, "_param.RData"))
SNP_XPASS = fit_bbj$SNP

#load Polypred results: SBAYESR EUR, SBAYESR target, Polyfun EUR

method = "Polypred"
result = fread(paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/polypred/",eth[1],"/",trait[l],"/SBayesR.snpRes"))
SNP_SBAYES_EUR = result$Name
result = fread(paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/polypred/",eth[i],"/",trait[l],"/SBayesR.snpRes"))
SNP_SBAYES_tar = result$Name
result = fread(paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/polypred/",eth[1],"/",trait[l],"/poly_fun"))
SNP_polyfun = result$SNP

#load prs-csx five ancestries results


file_out = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/XPASS/",eth[i],"/",trait[l],"/")
load(paste0(file_out, "_param.RData"))
SNP_XPASS = fit_bbj$SNP


#connect chr data into one
phi = c("1e+00","1e-02","1e-04","1e-06")

#for(i in 2:5){
out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/PRSCSx/",eth[i],"/",trait[l],"/")
#   for(l in 1:7){
v = 1
load(paste0(out.dir.prs,"sum_",eth[i],"_pst_eff_a1_b0.5_phi",phi[v],".rdata"))
snp.id = data$V2
load(paste0(out.dir.prs,"sum_EUR_pst_eff_a1_b0.5_phi",phi[v],".rdata"))
snp.id = data.frame(SNP = unique(c(snp.id,data$V2)))
snp.id = left_join(snp.id,snpinfo_mega_filter,by=c("SNP"="assay.name"))
BETA = data.frame(matrix(0,nrow=nrow(snp.id),ncol=4))
#load effect for target population
for(v in 1:4){
  load(paste0(out.dir.prs,"sum_",eth[i],"_pst_eff_a1_b0.5_phi",phi[v],".rdata"))
  temp.data = left_join(snp.id,data,by=c("SNP"="V2")) %>% 
    mutate(BETA = ifelse(is.na(V6),0,V6))
  BETA[,v] = temp.data %>% select(BETA)
}
BETA_tar = BETA
#load effect for EUR population
for(v in 1:4){
  load(paste0(out.dir.prs,"sum_",eth[1],"_pst_eff_a1_b0.5_phi",phi[v],".rdata"))
  temp.data = left_join(snp.id,data,by=c("SNP"="V2")) %>% 
    mutate(BETA = ifelse(is.na(V6),0,V6))
  BETA[,v] = temp.data %>% select(BETA)
}
BETA_EUR = BETA
data = cbind(snp.id,BETA_tar,BETA_EUR)
prs.snp = data[,-1]
out.dir.organize.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/organize_prs/PRSCSx/",eth[i],"/",trait[l],"/")
write.table(prs.snp,file = paste0(out.dir.organize.prs,"prs.file"),row.names = F,col.names = F,quote=F)
#   }



temp = fread( paste0(out.dir.organize.prs,"prs.file"))
# }

