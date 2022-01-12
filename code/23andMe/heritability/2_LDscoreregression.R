#calculate the LD-score for each ethnic using 1000 genome data
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])

eth = c("EUR","AFR","AMR","EAS","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")
data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"
library(data.table)
library(dplyr)
# for(i in 1:5){
#   for(l in 1:5){
    sid<-Sys.getenv('SLURM_JOB_ID')
    dir.create(paste0('/lscratch/',sid,'/test/'),showWarnings = FALSE)
    temp.dir = paste0('/lscratch/',sid,'/test/')
    if(l%in%c(3,4)){#continous traits
      sum.tar = as.data.frame(fread(paste0(data.dir,eth[i],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
      #load("/gpfs/gsfs11/users/zhangh24/multi_ethnic/data/cleaned/EUR/sumdat/morning_person_passQC_noNA_matchinfo.RData")
      sum.data = sum.tar %>% 
        mutate(MAF = ifelse(FREQ_A1<=0.5,FREQ_A1,1-FREQ_A1),
               N = N_control,
               se = SD,
               info =1,
               snpid  = rsid,
               Z = BETA/SD) %>% 
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
                    "/data/zhangh24/ldsc/ldsc.py ",
                    "--h2 ",temp.dir,"sum_align.sumstats.gz ",
                    "--ref-ld-chr ",ref.dir,eth[i],"_ldsc/ ",
                    "--w-ld-chr ",ref.dir,eth[i],"_ldsc/ ",
                    "--out ",temp.dir,"result "))
      result.folder = paste0("/data/zhangh24/multi_ethnic/result/cleaned/herit/",eth[i],"/",trait[l])
      system(paste0("less ",temp.dir,"result.log"))
      system(paste0("mv ",temp.dir,"result.log ",result.folder ))
      
    }else{
      sum.tar = as.data.frame(fread(paste0(data.dir,eth[i],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
      
      sum.data = sum.tar %>% 
        mutate(MAF = ifelse(FREQ_A1<=0.5,FREQ_A1,1-FREQ_A1),
               N = 1/(1/N_control+1/N_case),
               or = exp(BETA),
               se = SD,
               info =1,
               snpid  = rsid,
               Z = BETA/SD) %>% 
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
                    "/data/zhangh24/ldsc/ldsc.py ",
                    "--h2 ",temp.dir,"sum_align.sumstats.gz ",
                    "--ref-ld-chr ",ref.dir,eth[i],"_ldsc/ ",
                    "--w-ld-chr ",ref.dir,eth[i],"_ldsc/ ",
                    "--out ",temp.dir,"result "))
      result.folder = paste0("/data/zhangh24/multi_ethnic/result/cleaned/herit/",eth[i],"/",trait[l])
      system(paste0("mv ",temp.dir,"result.log ",result.folder ))
      
    }
    
#   }
# }
