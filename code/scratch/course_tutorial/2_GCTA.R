#Running GCTA to get the heritability
#Data were generated using 1000 Genomes European Data CHR 22 data
#Heritability was set up as 0.2
#Causal SNPs proportion: 5%

#Step 1, Compute the Genetic relationship matrix (GRM)
cur_dir <- "/data/BB_Bioinformatics/stat_gene_course/"
setwd(cur_dir)

system(paste0(cur_dir,"software/gcta64 ",
              "--bfile ",cur_dir,"data/chr22 ",
              "--make-grm ",
              "--out ",cur_dir,"result/chr22"))

#Step 2, Compute the Heritability
system(paste0(cur_dir,"software/gcta64 ",
              "--reml ",
              "--grm ",cur_dir,"result/chr22 ",
              "--pheno ",cur_dir,"result/phenotype.phen ",
              "--grm-cutoff 0.05 ",
              "--out ",cur_dir,"result/gcta_herit"))
