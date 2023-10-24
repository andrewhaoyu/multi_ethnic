#load the summary statistics

library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")

#load summary level statistics for overall
#Haoyupdate added the meta-analysis between iCOGs and OncoArray
standard_result <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF_Haoyuupdate.txt"))
colnames(standard_result)[1] <- "MarkerName"

bcac_result = fread("/data/zhangh24/ldsc/bcac_result.txt")
bcac_result = bcac_result[,-c(1:2)]
write.table(bcac_result,file = "/data/zhangh24/multi_ethnic/result/stat_gene_course/data/overall_bc",
            row.names = F, col.names = T, quote = F)
system(paste0("module load ldsc; ",
              "munge_sumstats.py ",
              "--sumstats /data/zhangh24/ldsc/bcac_result.txt ",
              "--out ",cur_dir,"result/ldsc_heritbcac ",
              "--merge-alleles /data/BB_Bioinformatics/software/ldsc_reference/eur_w_ld_chr/w_hm3.snplist ",
              " --chunksize 500000 ",
              "--signed-sumstats Z,0 --info-min 0.3 --maf-min 0.01"))

temp <- read.table(gzfile("/data/zhangh24/multi_ethnic/result/stat_gene_course/result/ldsc_heritbcac.sumstats.gz"))
