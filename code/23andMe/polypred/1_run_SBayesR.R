args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
v = as.numeric(args[[3]])
#j = as.numeric(args[[3]])
#m = as.numeric(args[[3]])
#i_rep = as.numeric(args[[4]])
#i1 = as.numeric(args[[5]])
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
sample_size <- as.data.frame(fread("/data/zhangh24/multi_ethnic/data/23_sample_size.csv",header=T))
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/1KGLD'),showWarnings = FALSE)
#copy the LD reference data to lscratch
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

#out.dir = paste0("/data/zhangh24/multi_ethnic/result/cleaned/clumping_result/",method,"/",eth[i],"/",trait[l],"/")

#prepare summary statistics

# input files
sum = as.data.frame(fread(paste0(data.dir,eth[i],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))

#match summary data with the reference data of SBayesR
LD_info = fread("/data/zhangh24/software/SBayesR/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_all_v3_50k.ldm.sparse.info", header = F)
colnames(LD_info)[1:6] = c("CHR", "rsid", "no1", "BP", "A1", "A2")

LD_info = LD_info %>% select(rsid) 

#combine summary data with SBayesR
#prepare the summary data format to be required format as SBayesR
summary_update = inner_join(sum, LD_info, by = "rsid")


bintrait <- c("any_cvd","depression",
              "iqb.sing_back_musical_note",
              "migraine_diagnosis",
              "morning_person")
contriat = c("heart_metabolic_disease_burden",
             "height")


if(trait[l]%in%bintrait){
  
  ma = summary_update %>% 
    select(rsid, A1, A2,FREQ_A1,  BETA, SD, P, N_control, N_case) %>% 
    rename(SNP = rsid,
           freq = FREQ_A1,
           b = BETA, 
           se = SD,
           p = P)
                    
}else{
  ma = summary_update %>% 
    mutate(N = N_control) %>% 
    select(rsid, A1, A2,FREQ_A1,  BETA, SD, P, N) %>% 
    rename(SNP = rsid,
           freq = FREQ_A1,
           b = BETA_align, 
           se = SD,
           p = P)
}

write.table(ma, file = paste0(temp.dir, "sumstats.ma"), row.names = F, col.names = T, quote = F)   


gctb_path = "/data/zhangh24/software/SBayesR/gctb_2.03beta_Linux/gctb"
#ref_path = paste0(LD.dir, "ukbEURu_hm3_chr",j,"_v3_50k.ldm.sparse")
summary_path = paste0(temp.dir, "sumstats.ma")
out_name =  paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/PRSCSx/",eth[i],"/",trait[l],"/SBayesR")

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

file_dir = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/PRSCSx/",eth[i],"/",trait[l])
files = dir(file_dir, pattern = paste0("SBayesR.snpRes"),full.names = T)
out_file = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/PRSCSx/",eth[i],"/",trait[l],"/SBayesR.snpRes")
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
  files = dir(file_dir, pattern = paste0("SBayesR.snpRes"),full.names = T)
  if(gamma_last<=1E-06){
    break
  }
}




# out.file = paste0("/data/zhangh24/multi_ethnic/result/cleaned/organize_prs/PRSCSx/")
# for(i in 1:5){
#   system(paste0("mkdir ",out.file,eth[i],"/"))
#   for(l in 1:7){
#     system(paste0("mkdir ",out.file,eth[i],"/",trait[l]))  
#   }
# }
