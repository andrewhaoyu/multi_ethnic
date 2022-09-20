args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
#v = as.numeric(args[[3]])
#j = as.numeric(args[[3]])
#m = as.numeric(args[[3]])
#i_rep = as.numeric(args[[4]])
#i1 = as.numeric(args[[5]])


eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")

library(data.table)
library(tidyverse)
sample_size <- as.data.frame(fread("/data/zhangh24/multi_ethnic/data/23_sample_size.csv",header=T))
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"
kg.dir = "/data/zhangh24/KGref_MEGA/GRCh37/"
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
system(paste0("cp ",kg.dir,eth[1],"/all_chr.bed ",temp.dir,eth[1],"all_chr.bed"))
system(paste0("cp ",kg.dir,eth[1],"/all_chr.bim ",temp.dir,eth[1],"all_chr.bim"))
system(paste0("cp ",kg.dir,eth[1],"/all_chr.fam ",temp.dir,eth[1],"all_chr.fam"))
system(paste0("cp ",kg.dir,eth[i],"/all_chr.bed ",temp.dir,eth[i],"all_chr.bed"))
system(paste0("cp ",kg.dir,eth[i],"/all_chr.bim ",temp.dir,eth[i],"all_chr.bim"))
system(paste0("cp ",kg.dir,eth[i],"/all_chr.fam ",temp.dir,eth[i],"all_chr.fam"))

#####
# Load library
library(XPASS)
library(data.table)
library(RhpcBLASctl)
#blas_set_num_threads(8)

#####
# parse the command line argument
#eval(parse(text=paste(commandArgs(trailingOnly = TRUE), collapse=";")))

### if not using command line, then please specify the parameters below.
# race = "AFR" # target population
# size = 1 # size for target population
# GA = 1
# rho = 1
# rep = 1
#####
# input files
summary_EUR = as.data.frame(fread(paste0(data.dir,eth[1],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))

#prepare the EUR data for XPASS format
summary_EUR_XPASS = summary_EUR %>% 
  select(SNP, N, Z, effect_allele, non_effect_allele) %>% 
  rename(A1 = effect_allele,
         A2 = non_effect_allele)
write.table(summary_EUR_XPASS, file = paste0(temp.dir,eth[1],"summary_eur"),
            row.names = F, col.names = T, quote = T)
#summary_EUR = paste0("/dcs04/nilanjan/data/wlu/XPASS/data/sumdata/EUR/sumdata-rho",rho,'-size4-rep',rep,'-GA',GA,'.txt') # auxilliary
summary_target = as.data.frame(fread(paste0(out.dir.sum,eth[i],"/summary_mega_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1))) 
summary_target_XPASS = summary_target %>% 
  select(SNP, N, Z, effect_allele, non_effect_allele) %>% 
  rename(A1 = effect_allele,
         A2 = non_effect_allele)
write.table(summary_target_XPASS, file = paste0(temp.dir,eth[i],"summary_tar"),
            row.names = F, col.names = T, quote = T)
path_to_summary_eur = paste0(temp.dir,eth[1],"summary_eur")
path_to_summary_tar = paste0(temp.dir,eth[i],"summary_tar")
#summary_target = paste0("/dcs04/nilanjan/data/wlu/XPASS/data/sumdata/",race,'/sumdata-rho',rho,'-size',size,'-rep',rep,'-GA',GA,'.txt') # target
# auxilliary
ref_gene_EUR = paste0(temp.dir,eth[1],"clump_ref_all_chr")
# target
ref_gene_target = paste0(temp.dir,eth[i],"clump_ref_all_chr")
file_out = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/xpass/rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)

#####
# XPASS
fit_bbj <- XPASS(file_z1 = path_to_summary_tar, file_z2 = path_to_summary_eur,
                 file_ref1 = ref_gene_target, file_ref2 = ref_gene_EUR,
                 # file_predGeno = ref_gene_pred, compPRS = T,
                 pop = "EUR", sd_method="LD_block", compPosMean = T,
                 file_out = file_out)
save(fit_bbj, file=paste0(file_out, "_param.RData"))
print("XPASS finished run")

system(paste0("cp ",cur.dir,eth[i],"/all_chr_test.mega.* ",temp.dir))

# predict
ref_gene_pred = paste0(temp.dir,"/all_chr_test.mega")
# output file prefix
load(paste0(file_out, "_param.RData"))
mu = fit_bbj$mu
# calculate PRS using plink (optional, comment the following out if not needed)
mu = mu[, c("SNP", "A1", "mu1", "mu2", "mu_XPASS1", "mu_XPASS2")]
write.table(mu,file = paste0(temp.dir,"prs_prep"),col.names = T,row.names = F,quote=F)
res = system(paste0("/data/zhangh24/software/plink2_alpha ",
                    "--score-col-nums 3,4,5,6 --threads 2 ",
                    "--score ",temp.dir,"prs_prep cols=+scoresums,-scoreavgs header no-mean-imputation ",
                    "--bfile ",ref_gene_pred,
                    " --out ",file_out,"_PRS"))





