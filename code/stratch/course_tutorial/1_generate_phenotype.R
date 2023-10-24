#select causal SNPs
library(data.table)
cur_dir <- "/data/BB_Bioinformatics/stat_gene_course/"
setwd(cur_dir)
snp <- fread("./data/chr22.bim")

causal_snp_idx <- sample(c(1:nrow(snp)),
                     0.05*nrow(snp),
                     replace = F)
causal_snp <- snp[causal_snp_idx,]
effect <- rnorm(nrow(causal_snp),0,sqrt(0.20/nrow(causal_snp)))
causal_snp_list <- data.frame(causal_snp$V2,effect)
write.table(causal_snp_list,
            file = paste0("./result/causal_snp_list"),
            row.names = F,
            col.names = F,
            quote = F)

system(paste0(cur_dir,"software/gcta64 ",
              "--bfile ",cur_dir,"data/chr22 ",
              "--simu-qt --simu-causal-loci ",cur_dir,"result/causal_snp_list ",
              "--simu-hsq 0.2 --simu-rep 1 --out ",cur_dir,"result/phenotype"))


temp <- fread(paste0(cur_dir,"result/phenotype.phen"))
var(temp$V3)

