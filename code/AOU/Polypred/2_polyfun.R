args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = 1
l = as.numeric(args[[1]])
#v = as.numeric(args[[3]])
#j = as.numeric(args[[3]])
#m = as.numeric(args[[3]])
#i_rep = as.numeric(args[[4]])
#i1 = as.numeric(args[[5]])
#"source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh; ",
#"mamba activate polyfun; ",
library(data.table, lib = "/home/zhangh24/R/4.1/library/")
library(dplyr, lib = "/home/zhangh24/R/4.1/library/") 
library(foreach, lib = "/home/zhangh24/R/4.1/library/")
library(iterators, lib = "/home/zhangh24/R/4.1/library/")
library(doParallel, lib = "/home/zhangh24/R/4.1/library/")
library(RcppZiggurat, lib = "/home/zhangh24/R/4.1/library/")
library(Rfast, lib = "/home/zhangh24/R/4.1/library/")



eth <- c("EUR","AFR","AMR")
trait <- c("height","bmi")
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
kg.dir = "/data/zhangh24/KGref_MEGA/GRCh37/"
system(paste0("cp ",kg.dir,eth[1],"/all_chr.bed ",temp.dir,eth[1],"all_chr.bed"))
system(paste0("cp ",kg.dir,eth[1],"/all_chr.bim ",temp.dir,eth[1],"all_chr.bim"))
system(paste0("cp ",kg.dir,eth[1],"/all_chr.fam ",temp.dir,eth[1],"all_chr.fam"))
data.dir = "/data/zhangh24/multi_ethnic/data/AOU_cleaned/"
#prepare summary statistics
#use the same preprocessing way as SBayesR
setwd("/data/zhangh24/multi_ethnic/")
# input files
#load sum data
sum = as.data.frame(fread(paste0(data.dir,eth[i],"/",trait[l],"_update.txt"),header=T))
colnames(sum)[1] = "rsid"
#match summary data with the reference data of SBayesR
LD_info = fread("/data/zhangh24/software/SBayesR/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_all_v3_50k.ldm.sparse.info", header = F)
colnames(LD_info)[1:6] = c("CHR", "rsid", "no1", "BP", "A1", "A2")

LD_info = LD_info %>% select(rsid) 
summary_update = inner_join(sum, LD_info, by = "rsid")
ma = summary_update %>% 
  rename(SNP = rsid,
         freq = A1_FREQ_allofus,
         b = BETA, 
         se = SE,
         p = P,
         BP  = pos37) %>% 
  mutate(Z = b/se,
         MAF = ifelse(freq<=0.5,freq, 1- freq)) %>% 
  select(SNP, CHR, BP, A1, A2,Z, N, MAF) 



write.table(ma, paste0(temp.dir, "sumstats"), row.names = F, col.names = T, quote = F)


############step1: munge the summary statistics to parquet format ###########
summary_path = paste0(temp.dir, "sumstats")

sum_file_align =  paste0(temp.dir,"sumstats_align")
system(
  paste0(
    "cd /data/zhangh24/software/polyfun;",
    "python munge_polyfun_sumstats.py ",
    "--sumstats ",summary_path," ",
    "--out ",sum_file_align," ",
    "--min-info 0 ",
    "--min-maf 0 ",
    "--keep-hla"
  )
)

#####################################################################
ref_gene_target = paste0(temp.dir,eth[i],"all_chr")

summary_path = paste0(temp.dir, "sumstats_align")

##################step 2:extract snp var from the database to the summary stats##############
system(
  paste0(
    "cd /data/zhangh24/software/polyfun;",
    "python extract_snpvar.py ",
    "--sumstats ",summary_path," ",
    "--out ",temp.dir,"sumstats_align_snpvar ",
    "--allow-missing ")
)
#####################################################################

#################step 3:Prepare the job file for polyfun##############


sum_file_align_snpvar = paste0(temp.dir,"sumstats_align_snpvar")

out_prefix_temp = paste0(temp.dir,"poly_fun_",trait[l])
N = median(ma$N)
system(paste0( "cd /data/zhangh24/software/polyfun; ",
               "python create_finemapper_jobs.py ",
               "--method susie ",
               "--sumstats ",sum_file_align_snpvar," ",
               "--n ",N," ",
               #"--non-funct ",
               "--geno ",ref_gene_target," ",
               #"--chr ",22, " "," ",
               "--memory 1"," ",
               "--max-num-causal 10 ",
               "--allow-missing ",
               "--out-prefix ",out_prefix_temp," ",
               "--jobs-file ",temp.dir,"job.sh "))
#####################################################################


#run the polyfun file in parallel
no.cores <- 8
#count the number of files


num <- as.integer(
  gsub(
    paste0(" ",temp.dir,"job.sh"),
    "", 
    system(paste0("wc -l ",temp.dir,"job.sh"),intern = T)
  )
)

#system(paste0("more ",temp.dir,"job.sh"))
# system(paste0("rm ",temp.dir,"job_split*"))
# system(paste0("ls ",temp.dir))
# system(paste0("wc -l ",temp.dir,"job_split00"))
# system(paste0("wc -l ",temp.dir,"job_split01"))
# system(paste0("wc -l ",temp.dir,"job.sh"))


system(paste0(
  "cd ",temp.dir,"; ",
  "split -l ",ceiling(num/no.cores), " -d job.sh job_split"))
#split the job file into half
system(paste0(
  "cd ",temp.dir,"; ",
  "split -l ",as.integer(num/no.cores), " -d job.sh job_split"))


registerDoParallel(no.cores)
result.list <- foreach(job.i = 1:no.cores)%dopar%{
  ID = paste0("0",c(0:(no.cores-1)))
  system(paste0("cd ",temp.dir,"; ",
                "sh ","job_split",ID[job.i]))
}
stopImplicitCluster()



#############step four: aggregate all the polyfun results into one
out_prefix = paste0("/data/zhangh24/multi_ethnic/result/AOU/prs/polypred/",eth[i],"/",trait[l],"/poly_fun")
system(
  paste0(
    "cd /data/zhangh24/software/polyfun; ",
    "python aggregate_finemapper_results.py ",
    "--out-prefix ",out_prefix_temp," ",
    "--sumstats ",sum_file_align_snpvar," ",
    "--out ",out_prefix," ", 
    "--adjust-beta-freq ",
    "--allow-missing-jobs "
  )
)


