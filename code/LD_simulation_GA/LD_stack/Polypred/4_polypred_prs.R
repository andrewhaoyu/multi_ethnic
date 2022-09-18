#get the prspolypred
## clear workspace
rm(list = ls())
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])
i1 =as.numeric(args[[5]])
library(data.table)
library(tidyverse)
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
system(paste0("cp ",cur.dir,eth[i],"/all_chr_test.mega.* ",temp.dir))



ref_gene_pred = paste0(temp.dir,"/all_chr_test.mega")
#calculate the PRS for SBayesR based on EUR 
file_dir = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[1],"/polypred")
files = dir(file_dir, pattern = ".snpRes", full.names = T)
file_out = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[1],"/polypred/rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)

file_coef = paste0(file_out,".snpRes")
#if the job converges, then nrow(file_coef) > 0
if((file_coef%in%files==T)&nrow(file_coef)>0){
  #load result
  sbayes_result = fread(paste0(file_out,".snpRes"),header = T)
  sbayes_result = sbayes_result %>% 
    mutate(chr.pos = paste0(Chrom,":",Position))


#load the original summary statistics to find 1KG ID
  summary <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_mega_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))
summary = summary %>% 
  mutate(chr.pos = paste0(CHR, ":", BP)) %>% 
  select(SNP,chr.pos)

#match sbayes R ID to 1KG ID
sbayes_result = left_join(sbayes_result, summary, by = "chr.pos")

assoc = sbayes_result %>% 
  rename(BETA = A1Effect) %>% 
  select(SNP, A1, BETA)
write.table(assoc,file = paste0(temp.dir,"prs_prep"),col.names = T,row.names = F,quote=F)
#calculate the PRS for SBayesR 
res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                    "--score-col-nums 3 --threads 2 ",
                    "--score ",temp.dir,"prs_prep cols=+scoresums,-scoreavgs header no-mean-imputation ",
                    "--bfile ",ref_gene_pred,
                    " --out ",file_out,"_PRS_SBayesR_EUR"))
}





#calculate the PRS for SBayesR based on target population


file_dir = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred")
files = dir(file_dir, pattern = ".snpRes", full.names = T)
file_out = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred/rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)

file_coef = paste0(file_out,".snpRes")
#load result
if((file_coef%in%files==T)&nrow(file_coef)>0){
sbayes_result = fread(paste0(file_out,".snpRes"),header = T)
sbayes_result = sbayes_result %>% 
  mutate(chr.pos = paste0(Chrom,":",Position))


#match sbayes R ID to 1KG ID
sbayes_result = left_join(sbayes_result, summary, by = "chr.pos")

assoc = sbayes_result %>% 
  rename(BETA = A1Effect) %>% 
  select(SNP, A1, BETA)
write.table(assoc,file = paste0(temp.dir,"prs_prep"),col.names = T,row.names = F,quote=F)
#calculate the PRS for SBayesR 
res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                    "--score-col-nums 3 --threads 2 ",
                    "--score ",temp.dir,"prs_prep cols=+scoresums,-scoreavgs header no-mean-imputation ",
                    "--bfile ",ref_gene_pred,
                    " --out ",file_out,"_PRS_SBayesR_tar"))
}

#calculate the PRS for Polyfun
#if(i ==1){
poly_fun_file_out = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[1],"/polypred/rho_",l,"_size_",4,"_rep_",i_rep,"_GA_",i1)
  polyfun_result = fread(poly_fun_file_out)  
  colnames(polyfun_result)[3] = "rsid"
  polyfun_result = polyfun_result %>% 
    mutate(chr.pos = paste0(CHR, ":", BP)) 
  polyfun_result = left_join(polyfun_result, summary, by = "chr.pos")
  
  assoc = polyfun_result %>% 
    rename(BETA = BETA_MEAN) %>% 
    select(SNP, A1, BETA)
  write.table(assoc,file = paste0(temp.dir,"prs_prep"),col.names = T,row.names = F,quote=F)
  #calculate the PRS for SBayesR
  res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                      "--score-col-nums 3 --threads 2 ",
                      "--score ",temp.dir,"prs_prep cols=+scoresums,-scoreavgs header no-mean-imputation ",
                      "--bfile ",ref_gene_pred,
                      " --out ",file_out,"_PRS_polyfun"))
#}

