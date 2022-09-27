args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
#j = as.numeric(args[[3]])
#m = as.numeric(args[[3]])
#i_rep = as.numeric(args[[4]])
#i1 = as.numeric(args[[5]])
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)
eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("HDL","LDL",
           "logTG",
           "TC")

sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
#copy the LD reference data to lscratch
system(paste0("mkdir ",temp.dir,"LD/"))
LD.dir = paste0(temp.dir,"LD/")
system(paste0("cp /data/zhangh24/software/SBayesR/ukbEURu_hm3_shrunk_sparse/*.bin ",LD.dir))
system(paste0("cp /data/zhangh24/software/SBayesR/ukbEURu_hm3_shrunk_sparse/*.info ",LD.dir))
data.dir = "/data/zhangh24/multi_ethnic/data/GLGC_cleaned/"


file_name  = rep("c", 22)
for(j in 1:22){
  file_name[j] = paste0(LD.dir, "ukbEURu_hm3_chr",j,"_v3_50k.ldm.sparse")
  
  
}
write.table(file_name, file = paste0(temp.dir, "ukbEURu_hm3_sparse_mldm_list.txt"),
            row.names = F, col.names = F, quote = F)
mldm_list_filename = paste0(temp.dir, "ukbEURu_hm3_sparse_mldm_list.txt")

#load sum data
sum = as.data.frame(fread(paste0(data.dir,eth[i],"/",trait[l],".txt"),header=T))
colnames(sum)[1] = "rsid"
#match summary data with the reference data of SBayesR
LD_info = fread("/data/zhangh24/software/SBayesR/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_all_v3_50k.ldm.sparse.info", header = F)
colnames(LD_info)[1:6] = c("CHR", "rsid", "no1", "BP", "A1", "A2")

LD_info = LD_info %>% select(rsid) 

#combine summary data with SBayesR
#prepare the summary data format to be required format as SBayesR
summary_update = inner_join(sum, LD_info, by = "rsid")

ma = summary_update %>% 
  rename(SNP = rsid,
         freq = A1_FREQ_1000G,
         b = BETA, 
         se = SE,
         p = P) %>% 
  select(SNP, A1, A2,freq,  b, se, p, N) 



write.table(ma, file = paste0(temp.dir, "sumstats.ma"), row.names = F, col.names = T, quote = F)   


gctb_path = "/data/zhangh24/software/SBayesR/gctb_2.03beta_Linux/gctb"
#ref_path = paste0(LD.dir, "ukbEURu_hm3_chr",j,"_v3_50k.ldm.sparse")
summary_path = paste0(temp.dir, "sumstats.ma")
out_name =  paste0("/data/zhangh24/multi_ethnic/result/GLGC/prs/polypred/",eth[i],"/",trait[l],"/SBayesR")

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

file_dir = paste0("/data/zhangh24/multi_ethnic/result/GLGC/prs/polypred/",eth[i],"/",trait[l])
files = dir(file_dir, pattern = paste0("SBayesR.snpRes"),full.names = T)
out_file = paste0("/data/zhangh24/multi_ethnic/result/GLGC/prs/polypred/",eth[i],"/",trait[l],"/SBayesR.snpRes")
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
    #model is not converged
    #out_file = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/Polypred/",eth[i],"/",trait[l],"/SBayesR.snpRes")
    result = NULL
    write.table(result, file = out_file,row.names = F, col.names = F, )
    break
  }
}
