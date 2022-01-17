args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)

library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)

eth <- c("EUR","AFR","AMR","EAS","SAS")
eth_group = c("european","african_american",
              "latino","east_asian","south_asian")

trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"
total.prs = length(pthres)
#load snp_info file to get corresponding im.id
load(paste0(data.dir,"snpinfo/snpinfo_mega.RData"))
snpinfo_mega_filter = snpinfo_mega %>% 
  filter(!is.na(im.data.id)) %>% 
  select(im.data.id,assay.name)
method = "TDLD_EBalleur"
#genearte index file to match columns in organize prs and parameters
# total = (length(eth)-1)*length(trait)*length(pthres)^2*length(r2_vec)*length(wc_base_vec)
# eth_vec = rep("c",total)
# trait_vec = rep("c",total)
# ptar_vec = rep(0,total)
# peur_vec = rep(0,total)
# col_vec = rep(0,total)
# r_ind_vec = rep(0,total)
# w_ind_vec = rep(0,total)

# temp = 1
# method = "TDLD"
# for(i in 2:length(eth)){
#   for(l in 1:length(trait)){
#     col_num = 1
#     for(r_ind in 1:length(r2_vec)){
#       for(w_ind in 1:length(wc_base_vec)){
#     for(k1 in 1:length(pthres)){
#       for(k2 in 1:length(pthres)){
#         col_vec[temp] = col_num
#         ptar_vec[temp] = k1
#         peur_vec[temp] = k2
#         trait_vec[temp] = trait[l]
#         eth_vec[temp] = eth[i]
#         r_ind_vec[temp] = r_ind
#         w_ind_vec[temp] = w_ind
#         col_num = col_num + 1
#         temp = temp+1
#       }
#     }
#       }
#     }
#    
#     }
#   }
# 
# index.match = data.frame(eth_vec,trait_vec,ptar_vec,peur_vec,col_vec,r_ind_vec,
#                          w_ind_vec)
# save(index.match,file = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/",method,"/index.match.rdata"))

# for(i in 1:length(eth)){
#   for(l in 1:length(trait)){
col_num = length(pthres)^2
#run through all the prs to get unique number of snps
rs.id.list = list()
out.dir.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/",method,"/",eth[i],"/",trait[l],"/")
for(i_post in 2:5){
  for(r_ind in 1:length(r2_vec)){
    for(w_ind in 1:length(wc_base_vec)){ 
      #just need the last p-thresholds to know the number of SNPs in particular r2 and window
      for(k1 in length(pthres):length(pthres)){
        for(k2 in length(pthres):length(pthres)){
          print(col_num)
          prs.file <- fread(paste0(out.dir.prs,col_num,"_",method,"_rind_",r_ind,"_wcind_",w_ind,"_ptar_",k1,"_peur_",k2),header=T)
          rs.id.list[[col_num]] = prs.file[,1,drop=F]
          col_num = col_num+length(pthres)^2
        }
      }
    }
  }
}

rs.id = unique(rbindlist(rs.id.list))
prs.snp = left_join(rs.id,snpinfo_mega_filter,by=c("SNP"="assay.name")) %>% 
  arrange(im.data.id)
rm(rs.id.list)
gc()

col_num = 1

total = 4*length(r2_vec)*length(wc_base_vec)*length(pthres)*length(pthres)
beta_mat = data.frame(matrix(nrow= nrow(prs.snp),ncol = total))

out.dir.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/",method,"/",eth[i],"/",trait[l],"/")

#generate prs following required format
col_num = 1
#column index to represent the current column position
col_ind = 1
#run through all the prs to get unique number of snps
out.dir.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/",method,"/",eth[i],"/",trait[l],"/")
library(caret)
for(i_post in 2:5){
  for(r_ind in 1:length(r2_vec)){
    for(w_ind in 1:length(wc_base_vec)){ 
      #filter out all the highly correlated columns
      beta_mat_temp = as.data.frame(matrix(nrow = nrow(prs.snp),ncol = length(pthres)^2))
      temp = 1
      for(k1 in 1:length(pthres)){
        for(k2 in 1:length(pthres)){
          print(col_num)
          prs.file <- fread(paste0(out.dir.prs,col_num,"_",method,"_rind_",r_ind,"_wcind_",w_ind,"_ptar_",k1,"_peur_",k2),header=T)
          prs.snp.temp = left_join(prs.snp,prs.file,by="SNP") %>% 
            mutate(BETA = ifelse(is.na(BETA),0,BETA)) %>%
            mutate(BETA = round(BETA,4)) 
          beta_mat_temp[,temp] = prs.snp.temp$BETA
          #colnames(prs.snp)[col_num+2] = paste0("BETA",col_num)
          col_num = col_num+1
          temp = temp + 1
        }
      }
      
      mtx = cor(beta_mat_temp)
      drop = findCorrelation(mtx,cutoff=0.98)
      drop = names(beta_mat_temp)[drop]
      beta_mat_temp_new = beta_mat_temp %>% 
        select(-all_of(drop))
      beta_mat[,col_ind:ncol(beta_mat_temp_new)] = beta_mat_temp_new
      col_ind = col_ind+ncol(beta_mat_temp_new)
    }
  }
  
}
#save(beta_mat,file = "/data/zhangh24/multi_ethnic/result/beta_mat.rdata")
beta_mat = beta_mat[,1:(col_ind-ncol(beta_mat_temp_new))]

prs.snp = cbind(prs.snp,beta_mat)
prs.snp = prs.snp[,-1]

out.dir.organize.prs <- paste0("/data/zhangh24/multi_ethnic/result/cleaned/organize_prs/",method,"/",eth_group[i],"/",trait[l],"/")
write.table(prs.snp,file = paste0(out.dir.organize.prs,"prs.file"),row.names = F,col.names = F,quote=F)
