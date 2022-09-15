#A summary statistics file with the following columns: SNP, CHR, BP, A1, A2, Z (Z-score), and (optionally) SNPVAR (per-SNP heritability, 
#which is proportional to prior causal probability). A1 is the effect allele (i.e. the sign of Z should match the A1 allele). 
#If the SNPVAR column is missing, you cannot perform functionally-informed fine-mapping.
#install.packages("data.table", lib = "/home/zhangh24/R/4.1/library/")
#install.packages("rlang", lib = "/home/zhangh24/R/4.1/library/")
#install.packages("dplyr", lib = "/home/zhangh24/R/4.1/library/")
#install.packages("vctrs", lib = "/home/zhangh24/R/4.1/library/")
#install.packages("doParallel", lib = "/home/zhangh24/R/4.1/library/")
#install.packages("foreach", lib = "/home/zhangh24/R/4.1/library/")
#install.packages("iterators", lib = "/home/zhangh24/R/4.1/library/")
#set up the environment of polyfun
#"source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh; ",
#"mamba activate polyfun; ",
rm(list = ls())
args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
m = as.numeric(args[[3]])
i_rep = as.numeric(args[[4]])
i1 =as.numeric(args[[5]])
library(data.table, lib = "/home/zhangh24/R/4.1/library/")
library(dplyr, lib = "/home/zhangh24/R/4.1/library/") 
library(foreach, lib = "/home/zhangh24/R/4.1/library/")
library(iterators, lib = "/home/zhangh24/R/4.1/library/")
library(doParallel, lib = "/home/zhangh24/R/4.1/library/")


eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
out.dir.sum <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/"
out.dir <-  "/data/zhangh24/multi_ethnic/result/LD_simulation_GA/LD_stack/"
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.bed ",temp.dir,eth[i],"clump_ref_all_chr.bed"))
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.bim ",temp.dir,eth[i],"clump_ref_all_chr.bim"))
system(paste0("cp ",cur.dir,eth[i],"/clump_ref_all_chr.fam ",temp.dir,eth[i],"clump_ref_all_chr.fam"))
#prepare summary statistics
#use the same preprocessing way as SBayesR
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
#"CHR", "BP", "SNP", "A1", "A2", "MAF", "N", "Z", "SNPVAR"
summary_update = summary_update %>% 
    mutate(BETA_align = ifelse(effect_allele == A1, BETA, -BETA),
           FQR_align = ifelse(effect_allele == A1, FQR_effect_allele,1 - FQR_effect_allele),
           SE = BETA/Z,
           MAF = ifelse(FQR_effect_allele<=0.5,FQR_effect_allele, 1-FQR_effect_allele))
    
    ma = summary_update %>% 
    select(RSID, CHR, BP, A1, A2,Z,N,MAF)
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
ref_gene_target = paste0(temp.dir,eth[i],"clump_ref_all_chr")

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
out_prefix_temp = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred/rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)

sum_file_align_snpvar = paste0(temp.dir,"sumstats_align_snpvar")

out_prefix_temp = paste0(temp.dir,"rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)

system(paste0( "cd /data/zhangh24/software/polyfun; ",
               "python create_finemapper_jobs.py ",
               "--method susie ",
               "--sumstats ",sum_file_align_snpvar," ",
               "--n 100000 ",
               #"--non-funct ",
               "--geno ",ref_gene_target," ",
              #"--chr ",22, " "," ",
               "--memory 1"," ",
               "--max-num-causal 10 ",
               "--allow-missing ",
               "--out-prefix ",out_prefix_temp," ",
               "--jobs-file ",temp.dir,"job.sh "))
#####################################################################



#count the number of files


num <- as.integer(
    gsub(
        paste0(" ",temp.dir,"job.sh"),
        "", 
        system(paste0("wc -l ",temp.dir,"job.sh"),intern = T)
        )
)

# system(paste0("rm ",temp.dir,"job_split*"))
# system(paste0("ls ",temp.dir))
# system(paste0("wc -l ",temp.dir,"job_split00"))
# system(paste0("wc -l ",temp.dir,"job_split01"))
# system(paste0("wc -l ",temp.dir,"job.sh"))


system(paste0(
    "cd ",temp.dir,"; ",
    "split -l ",ceiling(num/2), " -d job.sh job_split"))
#split the job file into half
system(paste0(
    "cd ",temp.dir,"; ",
    "split -l ",as.integer(num/2), " -d job.sh job_split"))

#run the polyfun file in parallel
no.cores <- 2
inner.size <- 2
registerDoParallel(no.cores)
result.list <- foreach(job.i = 1:inner.size)%dopar%{
    ID = c("00","01")
system(paste0("cd ",temp.dir,"; ",
              "sh ","job_split",ID[job.i]))
}
stopImplicitCluster()



#############step four: aggregate all the polyfun results into one
out_prefix = paste0("/data/zhangh24/multi_ethnic/result/LD_simulation_GA/",eth[i],"/polypred/rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1)
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





# shell_cmd = """python create_finemapper_jobs.py \
#                 --method susie \
#                 --sumstats {sumstats} \
#                 --n {size} \
#                 --geno {geno} \
#                 --chr {chr_num} \
#                 --memory 1 \
#                 --max-num-causal 10 \
#                 --allow-missing \
#                 --python3 python \
#                 --out-prefix {out_prefix} \
#                 --jobs-file {job_file} \
#                 2>&1
#                 """

