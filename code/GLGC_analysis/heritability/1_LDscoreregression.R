#calculate the LD-score for each ethnic using 1000 genome data
args = commandArgs(trailingOnly = T)
i = 1
l = as.numeric(args[[1]])

eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("HDL","LDL",
           "logTG",
           "TC")
data.dir = "/data/zhangh24/multi_ethnic/data/GLGC_cleaned/"
library(data.table)
library(dplyr)
# for(i in 1:5){
#   for(l in 1:5){
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test/'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')

sum.tar = as.data.frame(fread(paste0(data.dir,"EUR/",trait[l],".txt"),header=T))
  #load("/gpfs/gsfs11/users/zhangh24/multi_ethnic/data/cleaned/EUR/sumdat/morning_person_passQC_noNA_matchinfo.RData")
  sum.data = sum.tar %>% 
    mutate(MAF = as.numeric(ifelse(A1_FREQ_1000G<=0.5,A1_FREQ_1000G,1-A1_FREQ_1000G)),
           info =1,
           N = as.numeric(N),
           CHR = as.numeric(CHR),
           BP = as.numeric(POS_b37),
           snpid  = rsID,
           Z = as.numeric(BETA/SE),
           P = as.numeric(P)) %>% 
    select(snpid,CHR,BP,A1,A2,Z,P,info,MAF,N)
  #select(snpid,CHR,BP,A1,A2,or,se,P,info,MAF,N)
  fwrite(sum.data,file = paste0(temp.dir,"sum_data"),row.names = F,col.names = T,quote = F,sep = " ")
  #system( paste0("more ",temp.dir,"sum_data"))
  ref.dir = paste0("/data/zhangh24/KGref_MEGA/GRCh37/",eth[i],"/")
  system(paste0("module load ldsc; ",
                "/data/zhangh24/ldsc/munge_sumstats.py",
                " --sumstats ",temp.dir,"sum_data",
                " --out ",temp.dir,"sum_align",
                " --merge-alleles ",ref.dir,"snp_list",
                " --signed-sumstats Z,0",
                " --info-min 0.3 --maf-min 0.05; "))
  system(paste0("module load ldsc; ",
                "munge_sumstats.py",
                " --sumstats ",temp.dir,"sum_data",
                " --out ",temp.dir,"sum_align",
                " --merge-alleles ",ref.dir,"snp_list",
                " --signed-sumstats Z,0",
                " --info-min 0.3 --maf-min 0.05; "))
  
  system(paste0("module load ldsc; ",
                "/data/zhangh24/ldsc/ldsc.py ",
                "--h2 ",temp.dir,"sum_align.sumstats.gz ",
                "--ref-ld-chr ",ref.dir,eth[i],"_ldsc/ ",
                "--w-ld-chr ",ref.dir,eth[i],"_ldsc/ ",
                "--out ",temp.dir,"result "))
  result.folder = paste0("/data/zhangh24/multi_ethnic/result/GLGC/herit/",eth[i],"/",trait[l])
  #system(paste0("less ",temp.dir,"result.log"))
  system(paste0("mv ",temp.dir,"result.log ",result.folder ))
  
#   }
# }
