#goal: read UKB summary stat
library(readr)

# Path to your .csv.bgz file
#file path to T2D: /data/BB_Bioinformatics/ProjectData/UKB_sumstats/sumstats/E4_DM2.gwas.imputed_v3.both_sexes.tsv.bgz
#file path to Triglycerides (TG): /data/BB_Bioinformatics/ProjectData/UKB_sumstats/sumstats/30870_raw.gwas.imputed_v3.both_sexes.varorder.tsv.bgz
#file path to Total Cholesterol: /data/BB_Bioinformatics/ProjectData/UKB_sumstats/sumstats/30690_raw.gwas.imputed_v3.female.varorder.tsv.bgz
bgz_file_path <- "/data/BB_Bioinformatics/ProjectData/UKB_sumstats/sumstats/30690_raw.gwas.imputed_v3.female.varorder.tsv.bgz"

# Use gzfile to read the .csv.bgz file directly
df <- read.csv(gzfile(bgz_file_path),sep = "\t")
idx <- which(df$variant=="2:27730940:T:C")
df[idx,]
