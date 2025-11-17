load("./data/results_individual_analysis_genome.Rdata")
data = results_individual_analysis_genome
# install.packages("qqman")
# install.packages("dplyr")
library(qqman)
idx <- which(data$MAF >=0.05)
data_sub = data[idx,]
data_sub = cbind(data$CHR, data$POS,data$pvalue)
data_sub$SNP = paste0("SNP",c(1:nrow(data_sub)))
colnames(data_sub) = c("SNP","CHR","BP","P")
png(filename = "mahattan_result.png", width = 10, height = 8,
    unit = "in",res = 300)
manhattan(data_sub)
dev.off()
png(filename = "qq_result.png", width = 10, height = 8,
    unit = "in",res = 300)
qq(data_sub$P)
dev.off()


