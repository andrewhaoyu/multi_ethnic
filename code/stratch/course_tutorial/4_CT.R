#Construct the PRS using clumping and thresholding
#Data were generated using 1000 Genomes European population as the reference data
#Summary statistics of CHR 22 were provided
#PLINK 1.9 will be used for the clumping
#PLINK 2.0 will be used for score calculation
#R2 between PRS and Y will be calculated using R


sum_data_file = "/data/BB_Bioinformatics/stat_gene_course/data/EUR_sum_data"
ref_data = "/data/BB_Bioinformatics/stat_gene_course/data/1kg_eur_22/chr_22"
out_file = "/data/BB_Bioinformatics/stat_gene_course/result/LD_clump"
#clumping with windowsize 500kb, clumping r2 0.01, max p-value 1.0
res = system(paste0("/data/BB_Bioinformatics/stat_gene_course/software/plink ",
"--bfile ",ref_data," ",
"--clump ",sum_data_file," ",
"--clump-p1 1 ",
"--clump-r2 0.1  ",
"--clump-kb 500 ",
"--out ",out_file))


#load the clump result
LD_clump = fread(paste0(out_file,".clumped"))[,3,drop=F]
#load the summary statistics
EUR_sum = fread(sum_data_file)

library(dplyr)
#match the LD clumping results with summary statistics
prs_prep = left_join(LD_clump,EUR_sum, by = "SNP")
head(prs_prep)
#SNPs are ranked from the smallest p-value to largest p-value
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
for(k in 1:length(pthres)){
  prs_coeff = prs_prep %>% 
    filter(P<=pthres[k]) %>% 
    select(SNP, A1, BETA) 
  
  write.table(prs_coeff, file = 
                paste0("/data/BB_Bioinformatics/stat_gene_course/result/prs_coeff_",k),
              row.names = F,
              col.names = T,
              quote = F)
  geno_file = "/data/BB_Bioinformatics/stat_gene_course/data/prs_genotype/chr22_test"
  prs_coeff_file = paste0("/data/BB_Bioinformatics/stat_gene_course/result/prs_coeff_",k)
  prs_out = paste0("/data/BB_Bioinformatics/stat_gene_course/result/prs_",k)
  res <- system(paste0("/data/BB_Bioinformatics/stat_gene_course/software/plink2 ",
                       "--score ",prs_coeff_file," cols=+scoresums,-scoreavgs header no-mean-imputation  ",
                       "--bfile ",geno_file," --out ",prs_out))
}


#Evaluate the performance of 
#We have 20,000 people for tuning and validation purpose
#ID:10,001-11,000 will be used for the tuning dataset: select best p-value thresholding cutoff
#ID:11,001-12,000 will be used for the validation dataset: report the final performance
#read the outcome
y_out = fread(  paste0("/data/BB_Bioinformatics/stat_gene_course/data/y_out"))
y_tun = y_out[1:10000,"y"]
y_vad = y_out[10001:20000,"y"]
#create a vector to same the performance
r2_vec_tun = rep(0,length(pthres))
for(k in 1:length(pthres)){
  prs = fread(paste0("/data/BB_Bioinformatics/stat_gene_course/result/prs_",k,".sscore"))
  prs_tun = prs$SCORE1_SUM[1:10000]
  model = lm(y_tun$y~prs_tun)
  r2_vec_tun[k] = summary(model)$r.squared
}
#find best performance on the tuning dataset
idx.max = which.max(r2_vec_tun)
#evaluate it on the validation 
prs = fread(paste0("/data/BB_Bioinformatics/stat_gene_course/result/prs_",idx.max,".sscore"))
prs_vad = prs$SCORE1_SUM[10001:20000]
model = lm(y_vad$y~prs_vad)
r2 = summary(model)$r.squared
print(r2)
