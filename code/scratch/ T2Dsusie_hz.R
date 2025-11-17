#####################################
##Main Analysis
#####################################
# source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh
# mamba activate polyfun
# R --version
# library(data.table, lib = "/home/zhangh24/R/4.1/library/")
# library(dplyr, lib = "/home/zhangh24/R/4.1/library/") 
# library(foreach, lib = "/home/zhangh24/R/4.1/library/")
# library(iterators, lib = "/home/zhangh24/R/4.1/library/")
# library(doParallel, lib = "/home/zhangh24/R/4.1/library/")
# library(RcppZiggurat, lib = "/home/zhangh24/R/4.1/library/")
# library(Rfast, lib = "/home/zhangh24/R/4.1/library/")
#install.packages("vctrs")
#install.packages("vroom")
#install.packages("tidyverse")
#install.packages("data.table")
#install.packages("aws.s3")
library(vctrs)
library(vroom)
library(data.table)
library(tidyverse)
library(aws.s3)
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')

#############################
#############################
trait_vec = c("T2D","TG","TC")
trait = trait_vec[3]
if(trait=="T2D"){
  file_path <- "/data/BB_Bioinformatics/DG/t2d/E4_DM2.gwas.imputed_v3.both_sexes.tsv.bgz"
}else if(trait == "TG"){
  trigl_path <- "/data/BB_Bioinformatics/ProjectData/UKB_sumstats/sumstats/30870_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz"
}else{
  chol_path <- "/data/BB_Bioinformatics/ProjectData/UKB_sumstats/sumstats/30690_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz" 
}
# 
#

#############################
#############################
summ <- vroom(chol_path)
sum <- summ %>% filter(startsWith(variant, "2:"))
dim(sum)
setDT(sum)
# Split the 'variant' column using data.table's tstrsplit function
sum[, c("chromosome", "position", "allele1", "allele2") := tstrsplit(variant, ":", fixed = TRUE)]
head(sum)
#############################
#############################
sum2 <- sum
sum2$chromosome <- as.numeric(sum2$chromosome)
sum2$position <- as.numeric(sum2$position)
sum2 <- mutate(sum2, A1 = ifelse(minor_allele == allele1, allele1, allele2))
sum2 <- mutate(sum2, A2 = ifelse(minor_allele == allele1, allele2, allele1))
head(sum2)
dim(sum2)
#############################
#############################

# bucket <- "broad-alkesgroup-ukbb-ld"
# file_path <- "UKBB_LD/chr2_25000001_28000001.gz"
# aws.s3::save_object(file_path, bucket = bucket, file ="/data/BB_Bioinformatics/DG/t2d/ld.gz")
# list.files("/data/BB_Bioinformatics/DG")
# ld <- vroom("/data/BB_Bioinformatics/DG/t2d/ld.gz")
# head(ld)

# file_path <- "UKBB_LD/chr2_25000001_28000001.npz"
# aws.s3::save_object(file_path, bucket = bucket, file ="/data/BB_Bioinformatics/DG/t2d/ld.npz")

# joined <- inner_join(sum2, ld[c("chromosome", "position", "rsid", "allele1", "allele2")], by = c("chromosome", "position", "allele1", "allele2"))
# joined <- joined %>% select(variant, minor_AF, n_complete_samples, beta, se, tstat, pval, chromosome, position, rsid, A1, A2)

#load the UKBB EUR data
ukb_bim <- as.data.frame(fread("/data/BB_Bioinformatics/ProjectData/UKBB_sub/genotype/all_data/EUR/all_chr.bim"))
colnames(ukb_bim) <- c(c("chromosome", "rsid", "nothing","position",  "allele1", "allele2"))


joined <- inner_join(sum2, ukb_bim[c("chromosome", "position", "rsid")], by = c("chromosome", "position"))

joined <- joined %>% 
  rename(SNP = rsid,
         b = beta, 
         p = pval,
         CHR = chromosome,
         BP  = position, 
         N = n_complete_samples,
         MAF = minor_AF) %>% 
  mutate(Z = b/se) %>% 
  select(SNP, CHR, BP, A1, A2, Z, N, MAF) 

#temp.dir <- "/data/BB_Bioinformatics/DG/t2d/"
write.table(joined, paste0(temp.dir, "stats"), row.names = F, col.names = T, quote = F)
#############################
#############################

###########################################
##step1: munge the summary statistics to parquet format
###########################################
summary_path = paste0(temp.dir, "stats")
sum_file_align =  paste0(temp.dir,"sumstats_align")

system(
  paste0(
    "cd /data/BB_Bioinformatics/software/polyfun;",
    "python munge_polyfun_sumstats.py ",
    "--sumstats ",summary_path," ",
    "--out ",sum_file_align," ",
    "--min-info 0 ",
    "--min-maf 0 ",
    "--keep-hla"
  )
)

summary_path = paste0(temp.dir, "sumstats_align")
###########################################
##step 2:extract snp var from the database to the summary stat
###########################################
system(
  paste0(
    "cd /data/BB_Bioinformatics/software/polyfun;",
    "python extract_snpvar.py ",
    "--sumstats ",summary_path," ",
    "--out ",temp.dir,"sumstats_align_snpvar ",
    "--allow-missing ")
)

###########################################
##step 3:Prepare the job file for polyfun
###########################################
sum_file_align_snpvar = paste0(temp.dir,"sumstats_align_snpvar")
temp = fread(sum_file_align_snpvar)
out_prefix_temp = paste0(temp.dir,"poly_fun_",trait)



N <- 343992
ref_gene_target <- "/data/BB_Bioinformatics/ProjectData/UKBB_sub/genotype/all_data/EUR/all_chr"
system(paste0( "cd /data/BB_Bioinformatics/software/polyfun; ",
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
#system("more /lscratch/7950833/test/job.sh")
#system("head -1 /lscratch/7950833/test/job.sh > /lscratch/7950833/test/test.sh")
#system(paste0("head -1 ",temp.dir,"job.sh > ",temp.dir,"test.sh"))
#system(paste0("sh ",temp.dir,"test.sh"))
#system(paste0("more ",temp.dir,"test.sh"))
system(paste0("sh ",temp.dir,"job.sh"))

#############step four: aggregate all the polyfun results into one
#out_prefix_temp = paste0(temp.dir,"poly_fun_",trait)
#Devika: change the out_prefix to your output directory
out_prefix = paste0("/data/BB_Bioinformatics/polyfun_test/poly_fun_",trait)
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




# system(paste0("mamba activate polyfun-lock;",
#   "cd /data/BB_Bioinformatics/software/polyfun; ",
#               "python finemapper.py ",
#               "--ld ", "/data/BB_Bioinformatics/DG/t2d/ld ",
#               "--sumstats ", sum_file_align_snpvar, " ",
#               "--n ", N, " ",
#               "--chr ", "2 ",
#               "--start ", "27719900 ",
#               "--end ", "27746463 ",
#               "--method susie ",
#               "--susie-resvar 1 ",
#               "--max-num-causal 10 ", 
#               "--out ", "/data/BB_Bioinformatics/DG/t2d/nonfuncdiabetes.gz"))



# missingsnp_path <- "/data/BB_Bioinformatics/DG/sumstats_align_snpvar.miss.gz"
# missingsnp <- vroom(missingsnp_path)
# missingsnp
# 
# nonfunc_path <- "/data/BB_Bioinformatics/DG/t2d/nonfuncdiabetes.gz"
# nonfunc <- vroom(nonfunc_path)
# nonfunc
# 
# func_path <- "/data/BB_Bioinformatics/DG/t2d/funcdiabetes.gz"
# func <- vroom(func_path)
# head(func)
