#Use LDSC to estimate heritability
cur_dir <- "/data/zhangh24/multi_ethnic/result/stat_gene_course/"
#look at the data
#the data is based on breast cancer GWAS
bcac_overall <- fread(paste0(cur_dir,"data/overall_bc"))
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
              "--merge-alleles /data/BB_Bioinformatics/software/ldsc_example_data/w_hm3.snplist ",
              "--chunksize 500000 ",
              "--signed-sumstats Z,0 --info-min 0.3 --maf-min 0.01"))

munge_result = paste0(cur_dir,"result/ldsc_herit_overall.sumstats.gz")

system(paste0("module load ldsc; ",
       "ldsc.py ",
       "--h2 ", munge_result," ",
       "--ref-ld-chr ", cur_dir,"data/eur_w_ld_chr/ ",
       "--w-ld-chr ", cur_dir,"data/eur_w_ld_chr/ ",
       "--out ", cur_dir,"result/h2_overall "))
#observed scale heritablity is 0.4777


