#Goal: prepare summary statistics for mega SNPs
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
i1 = as.numeric(args[[2]])
library(data.table)
library(tidyverse)
cur.dir = "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir.sum = "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
eth = c("EUR", "AFR", "AMR", "EAS", "SAS")

setwd("/data/zhangh24/multi_ethnic/")
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")

snp.infor = snp.infor.match %>% 
  rename(SNP=id) %>% 
  select(SNP,rs_id)
mega.list <- as.data.frame(fread(paste0(cur.dir,"mega-hm3-rsid.txt"), header = F))
colnames(mega.list) = "rs_id"

#for(i in 1:5){
  #for(i1 in 1:5){
    for(l in 1:3){
      for(m in 1:4){
        for(i_rep in 1:10){
          summary = as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_out_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))) 
          summary = left_join(summary,snp.infor,by="SNP")
          #subset SNPs to mega SNPs
          summary_match = inner_join(summary,mega.list,by="rs_id")
          #A1 is the effect allele
          #A2 is the non-effect allele; generate A2 as non-effect-allele
          summary_match_update = summary_match %>% 
            separate(SNP,into = c("non_im","non_im2","first_allele","second_allele"),
                     remove = F) %>% 
            mutate(A2 =ifelse(A1==second_allele,first_allele,second_allele)) 
          
          summary_match_update = summary_match_update %>% 
            select(SNP, rs_id, CHR, BP, A1, A2, BETA, STAT, P, NMISS) %>% 
            rename(effect_allele = A1,
                   non_effect_allele = A2,
                   N = NMISS,
                   Z = STAT)
          write.table(summary_match_update, file = paste0(out.dir.sum,eth[i],"/summary_mega_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1),
                      row.names = F, col.names = T, quote = F)    
        }
      }
    }
  #}  
 # }
  
