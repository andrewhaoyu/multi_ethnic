

suppressMessages(library('plink2R'))

genos = read_plink(paste("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/1000G_final/chr2122"), impute="avg")
pc1000G <- prcomp(genos$bed, center=TRUE, scale=TRUE, rank=10)
saveRDS(pc1000G, "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/pc1000G_5.rds")

writeLines(rownames(pc1000G$rotation), "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/pc1000G_5.snps")
saveRDS(pc1000G$rotation, "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/pc1000G_5.loadings")
saveRDS(pc1000G$x, "/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/pc1000G_5.pcscores")

scov <- cov(t(genos$bed))
saveRDS(scov,"/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/hdpca/scov.rds")
ev <- eigen(scov)
saveRDS(ev,"/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/hdpca/ev.rds")


genomes_sd <- apply(genos$bed, MARGIN=2, sd)
genomes_mean <- apply(genos$bed, MARGIN=2, mean)
save(genomes_sd, genomes_mean, file="/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/hdpca/geno1000G.sdmean.RData")


library(hdpca)

#First estimate the number of spikes and then adjust test scores based on that
train.eval<-ev$values
n<-nrow(genos$fam)
p<-nrow(genos$bim)
m<-select.nspike(train.eval,p,n,n.spikes.max=50,evals.out=FALSE)$n.spikes ## 33
save(train.eval, n, p, m, file="/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/models/hdpca/parameters.RData")
