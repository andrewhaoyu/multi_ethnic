#Goal: preprocess the summary statistics to add the allele frequency to the summary stat
rm(list = ls())
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i1 =as.numeric(args[[4]])
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
#load 1kg snp list
setwd("/data/zhangh24/multi_ethnic/")
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
snp_infor_match_sub = snp.infor.match %>% 
  select(id, a0, a1, eth[i]) %>% 
  rename(FRQ_A1 = eth[i])
snp_infor_match_sub = snp_infor_match_sub %>% distinct()

for(i_rep in 1:10){
  # input files
  summary <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_mega_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))
  #remove duplicated rows
  summary = summary %>% distinct(SNP, .keep_all = TRUE)
  summary_update = left_join(summary, snp_infor_match_sub, by = c("SNP" = "id"))
  #create the allele frequency for effect allele
  summary_update = summary_update %>% 
    mutate(
      FQR_effect_allele = case_when(
        effect_allele==a1 ~ FRQ_A1,
        effect_allele==a0 ~ 1-FRQ_A1
      )
    )
  summary_update = summary_update %>% 
    select(c(colnames(summary),"FQR_effect_allele"))
  # write.table(summary_update, file =   paste0(out.dir.sum,eth[i],"/summary_mega_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1),
  #             row.names = F, col.names = T, quote = F)
  write.table(summary_update, file =   paste0(out.dir.sum,eth[i],"/test"),
              row.names = F, col.names = T, quote = F)
  
}
              
          
              
              
              
              
              
            