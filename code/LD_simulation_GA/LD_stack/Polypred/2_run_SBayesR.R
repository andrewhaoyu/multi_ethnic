# run_XPASS.R runs XPASS and calculate PRS scores (by both XPASS and plink,
# comment the last few lines out if plink is not needed). Note that the way
# to call plink might vary. (It's plink2 --... for me.)
#
# Run with:
# R CMD BATCH --vanill '--args race="AFR" size=1 GA=1 rho=2 rep=1' run_XPASS.R
#
# The section of "input files" controls the input and output file paths.
# Output files have the prefix paste0("/dcs04/nilanjan/data/wlu/XPASS/XPASS_result/", race, "/XPASS-rho",rho,'-size',size,'-rep',rep,'-GA',GA)

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

setwd("/data/zhangh24/multi_ethnic/")
# input files
summary <- as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_mega_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)))
summary = summary %>% 
  mutate(chr.pos = paste0(CHR, ":", BP))
#match summary data with the reference data of SBayesR
LD_info = fread("/data/zhangh24/software/SBayesR/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_all_v3_50k.ldm.sparse.info", header = F)
colnames(LD_info)[1:6] = c("CHR", "RSID", "no1", "BP", "A1", "A2")

LD_info = LD_info %>% select(CHR, RSID, BP, A1, A2) %>% 
  mutate(chr.pos = paste0(CHR, ":", BP)) %>% 
  select(RSID, chr.pos, A1, A2)

#combine summary data with SBayesR
#align the effect to be the same as the reference of SBayesR
#prepare the summary data format to be required format as SBayesR
summary_update = inner_join(summary, LD_info, by = "chr.pos")

summary_update = summary_update %>% 
  mutate(BETA_align = ifelse(effect_allele == A1, BETA, -BETA),
         FQR_align = ifelse(effect_allele == A1, FQR_effect_allele,1 - FQR_effect_allele),
         SE = BETA/Z)

#for(j in 1:22){
  #copy the reference data
  system(paste0("mkdir ",temp.dir,"LD/"))
  LD.dir = paste0(temp.dir,"LD/")
  system(paste0("cp /data/zhangh24/software/SBayesR/ukbEURu_hm3_shrunk_sparse/*.bin ",LD.dir))
  system(paste0("cp /data/zhangh24/software/SBayesR/ukbEURu_hm3_shrunk_sparse/*.info ",LD.dir))
  system(paste0("ls ",LD.dir))
  #create reference file list
  file_name  = rep("c", 22)
  for(j in 1:22){
    file_name[j] = paste0(LD.dir, "ukbEURu_hm3_chr",j,"_v3_50k.ldm.sparse")
    
    
  }
  write.table(file_name, file = paste0(temp.dir, "ukbEURu_hm3_sparse_mldm_list.txt"),
              row.names = F, col.names = F, quote = F)
  mldm_list_filename = paste0(temp.dir, "ukbEURu_hm3_sparse_mldm_list.txt")
  #prepare the summary data format to be required format as SBayesR
  ma = summary_update %>% 
    select(RSID, A1, A2,FQR_align,  BETA_align, SE, P, N) %>% 
    rename(SNP = RSID,
           freq = FQR_align,
           b = BETA_align, 
           se = SE,
           p = P)
  write.table(ma, file = paste0(temp.dir, "sumstats.ma"), row.names = F, col.names = T, quote = F)   
  
  
  
  gctb_path = "/data/zhangh24/software/SBayesR/gctb_2.03beta_Linux/gctb"
  #ref_path = paste0(LD.dir, "ukbEURu_hm3_chr",j,"_v3_50k.ldm.sparse")
  summary_path = paste0(temp.dir, "sumstats.ma")
  out_name =  paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred/rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)
  system(paste0(gctb_path," --sbayes R ",
                "--mldm ", mldm_list_filename," ",
                "--pi 0.95,0.02,0.02,0.01 ",
                "--gamma 0.0,0.01,0.1,1 ",
                "--gwas-summary ",summary_path," ",
                "--chain-length 10000 ",
                "--burn-in 4000 ",
                "--out-freq 10 ",
                "--thin 10 ",
                "--out ",out_name))
  #if SBayesR doesn't converge, change the gamma last to gamma_last/2 until convergence
  file_dir = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred")
  files = dir(file_dir, pattern = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".snpRes"),full.names = T)
  out_file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred/rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".snpRes")
  gamma_last = 1
  while(out_file%in%files==F){
    gamma_last = gamma_last/2
    system(paste0(gctb_path," --sbayes R ",
                  "--mldm ", mldm_list_filename," ",
                  "--pi 0.95,0.02,0.02,0.01 ",
                  "--gamma 0.0,0.01,0.1,",gamma_last," ",
                  "--gwas-summary ",summary_path," ",
                  "--chain-length 10000 ",
                  "--burn-in 4000 ",
                  "--out-freq 10 ",
                  "--thin 10 ",
                  "--out ",out_name))
    files = dir(file_dir, pattern = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1),full.names = T)
    if(gamma_last<=1E-06){
      result = NULL
      write.table(result, file = out_file,row.names = F, col.names = F, )
      break
    }
  }
#   
# #}     
#   mis_vec_list = list()   
#   
#   mis_temp =  1
#   
#   for(i in 1:1){
#     file_dir = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred")
#     
#     for(l in 1:3){
#       for(m in 4:4){
#         for(i_rep in 1:10){
#           for(i1 in 1:5){
#             files = dir(file_dir, pattern = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".snpRes"),full.names = T)      
#             out_file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred/rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".snpRes")
#             if(out_file%in%files==F){
#               mis_vec_list[[mis_temp]] = data.frame(i, l, m, i_rep, i1)
#               mis_temp = mis_temp + 1
#             }
#             
#           }
#         }
#       }
#     }
#   }
#   
#   for(i in 2:5){
#     file_dir = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred")
#     
#     for(l in 1:3){
#       for(m in 1:4){
#         for(i_rep in 1:10){
#           for(i1 in 1:5){
#             files = dir(file_dir, pattern = paste0("rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".snpRes"),full.names = T)      
#             out_file = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred/rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".snpRes")
#             if(out_file%in%files==F){
#               mis_vec_list[[mis_temp]] = data.frame(i, l, m, i_rep, i1)
#               mis_temp = mis_temp + 1
#             }
#             
#           }
#         }
#       }
#     }
#   }
#   
#   mis_vec = rbindlist(mis_vec_list)
#   run_job = rep("c", nrow(mis_vec))
#   for(k in 1:nrow(mis_vec)){
#     run_job[k] = paste0("Rscript /gpfs/gsfs11/users/zhangh24/multi_ethnic/code/LD_simulation_GA/LD_stack/Polypred/2_run_SBayesR.R ",
#                          mis_vec[k,1]," ",mis_vec[k,2]," ",mis_vec[k,3]," ",mis_vec[k,4]," ",mis_vec[k,5])
#   }
#   write.table(run_job, file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation_GA/LD_stack/Polypred/rerun_SBayesR.sh"),
#               row.names = F, col.names = F, quote = F)
#       