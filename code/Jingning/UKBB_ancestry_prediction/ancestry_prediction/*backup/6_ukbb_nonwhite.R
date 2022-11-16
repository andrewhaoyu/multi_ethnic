

suppressMessages(library('plink2R'))

genos = read_plink(paste("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_final/chr2122"), impute="avg")
pc1000G <- prcomp(genos$bed, center=TRUE, scale=TRUE, rank=10)
saveRDS(pc1000G, "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/pc1000G_5.rds")

writeLines(rownames(pc1000G$rotation), "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/pc1000G_5.snps")
saveRDS(pc1000G$rotation, "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/pc1000G_5.loadings")
saveRDS(pc1000G$x, "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/pc1000G_5.pcscores")
