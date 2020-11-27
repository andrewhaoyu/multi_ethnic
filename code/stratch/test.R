library(bigsnpr)
setwd("/data/zhangh24/multi_ethnic/result")
snp_readBed("tmp-data/public-data.bed")
obj.bigSNP <- snp_attach("tmp-data/public-data.rds")
str(obj.bigSNP, max.level = 2, strict.width = "cut")
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection - 1
NCORES <- nb_cores()
# Check some counts for the 10 first variants
big_counts(G, ind.col = 1:10)

sumstats <- bigreadr::fread2("tmp-data/public-data-sumstats.txt",
                             select = c(1:6, 10))
str(sumstats)

set.seed(1)
ind.train <- sample(nrow(G), 400)
ind.test <- setdiff(rows_along(G), ind.train)

names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "p")
map <- obj.bigSNP$map[,-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats, map)

beta <- info_snp$beta
lpval <- -log10(info_snp$p)

NCORES = 2
all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train,
                              lpS = lpval)
attr(all_keep, "grid")
