#Use LDSC to estimate heritability
#Data were obtained from the GWAS summary statistics of breast cancer
#Three different traits were included: overall breast cancer risk, Luminal A, Triple negative
#More background of the GWAS can be found in: https://www.nature.com/articles/s41588-020-0609-2
cur_dir <- "/data/BB_Bioinformatics/stat_gene_course/"
#look at the data
#the data is based on breast cancer GWAS
bcac_overall <- fread(paste0(cur_dir,"data/overall_bc"))
head(bcac_overall)
#snpid; if the SNP has rs id, then snpid is rsid, otherwise snpid is chr:position
#CHR: chromosome
#bp: GRCh37 (hg19)
#A1: effect_allele
#A2: non_effect_allele
#Z: Z-statistics
#P: P-value
#info: imputation quality score
#MAF: minor allele frequency
#N: effective-sample size; N was calculated as 1/(var(beta)*2*f*(1-f))
#effective sample size is needed if one wants to calculate the logit-scale genetic variance
system(paste0("module load ldsc; ",
              "munge_sumstats.py ",
              "--sumstats ",cur_dir,"data/overall_bc ",
              "--out ",cur_dir,"result/ldsc_herit_overall ",
              "--merge-alleles ",cur_dir,"data/eur_w_ld_chr/w_hm3.snplist ",
              "--chunksize 500000 ",
              "--signed-sumstats Z,0 --info-min 0.3 --maf-min 0.01"))

munge_result = paste0(cur_dir,"result/ldsc_herit_overall.sumstats.gz")

#the frality scale heritablity is 0.4777
system(paste0("module load ldsc; ",
       "ldsc.py ",
       "--h2 ", munge_result," ",
       "--ref-ld-chr ", cur_dir,"data/eur_w_ld_chr/ ",
       "--w-ld-chr ", cur_dir,"data/eur_w_ld_chr/ ",
       "--out ", cur_dir,"result/h2_overall "))



#genetic correlation calculation for luminal A and triple negative breast cancer subtypes
#munge the luminal A and triple negative summary statistics
lua <- fread(paste0(cur_dir,"data/lua_bc"))
head(lua)
system(paste0("module load ldsc; ",
              "munge_sumstats.py ",
              "--sumstats ",cur_dir,"data/lua_bc ",
              "--out ",cur_dir,"result/ldsc_herit_lua ",
              "--merge-alleles ",cur_dir,"data/eur_w_ld_chr/w_hm3.snplist ",
              "--chunksize 500000 ",
              "--signed-sumstats Z,0 --info-min 0.3 --maf-min 0.01"))

munge_result_lua = paste0(cur_dir,"result/ldsc_herit_lua.sumstats.gz")

#munge the luminal A and triple negative summary statistics
system(paste0("module load ldsc; ",
              "munge_sumstats.py ",
              "--sumstats ",cur_dir,"data/tn_bc ",
              "--out ",cur_dir,"result/ldsc_herit_tn ",
              "--merge-alleles ",cur_dir,"data/eur_w_ld_chr/w_hm3.snplist ",
              "--chunksize 500000 ",
              "--signed-sumstats Z,0 --info-min 0.3 --maf-min 0.01"))

munge_result_tn = paste0(cur_dir,"result/ldsc_herit_tn.sumstats.gz")

#calculate genetic correlation
system(paste0("module load ldsc; ",
              "ldsc.py ",
              "--rg ", munge_result_lua,",",munge_result_tn," ",
              "--ref-ld-chr ", cur_dir,"data/eur_w_ld_chr/ ",
              "--w-ld-chr ", cur_dir,"data/eur_w_ld_chr/ ",
              "--out ", cur_dir,"result/rg_lua_tn "))
#genetic correlation 0.4829 (s.e. 0.0512)


#stratified LD-score regression using baseline annotation
system(paste0("module load ldsc; ",
              "ldsc.py ",
              "--h2 ", munge_result," ",
              "--ref-ld-chr ", cur_dir,"data/1000G_Phase3_baselineLD_ldscores/baselineLD. ",
              "--w-ld-chr ", cur_dir,"data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. ",
              "--overlap-annot  ",
              "--frqfile-chr ", cur_dir,"data/1000G_Phase3_frq/1000G.EUR.QC. ",
              "--out ", cur_dir,"result/h2_sldsc "))
library(data.table)
enrichment_result = fread("/data/BB_Bioinformatics/stat_gene_course/result/h2_sldsc.results")
which.max(enrichment_result$Enrichment)
enrichment_result[13,]
