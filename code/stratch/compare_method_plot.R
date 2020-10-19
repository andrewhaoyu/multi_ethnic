setwd("/Users/zhangh24/Desktop/compare_method_plot")
load("data_for_plot.rdata")
source("theme_publication.R")
library(RColorBrewer)
colourCount = length(unique(LD.clump.result.plot$method_vec))
getPalette = colorRampPalette(brewer.pal(9, "Paired"))
method <- as.character(LD.clump.result.plot$method_vec)
idx <- which(method=="P+T")
method[idx] <- "Population-specific P+T"
idx <- which(method=="LDpred")
method[idx] <- "Population-specific LDpred"
idx <- which(method=="EUR P+T PRS")
method[idx] <- "EUR P+T"
idx <- which(method=="EUR LDpred PRS")
method[idx] <- "EUR LDpred"
idx <- which(method=="2DLD")
method[idx] <- "2DP+T"
idx <- which(method=="2DLD-EB")
method[idx] <- "2DP+T-EB"

method.fac <- factor(method,
                     levels = c("Population-specific P+T",
                                "Population-specific LDpred",
                                "EUR P+T",
                                "EUR LDpred",
                                "SNPs in EUR P+T + target coefficients",
                                "SNPs in EUR P+T + EB coefficients",
                                "2DP+T",
                                "2DP+T-EB",
                                "ME-Bayes"))
LD.clump.result.plot$method_vec <- method.fac


p <- ggplot(LD.clump.result.plot,aes(x= sample_size,y=r2.vec,group=method_vec))+
  geom_bar(aes(fill=method_vec),
           stat="identity",
           position = position_dodge())+
  #geom_point(aes(color=method_vec))+
  theme_Publication()+
  ylab("R2")+
  xlab("Sample size for minority group")+
  labs(fill = "Method")+
  facet_grid(vars(cau_vec),vars(eth.vec))+
  scale_fill_manual(values = getPalette(colourCount))+
  ggtitle("Prediction comparasion (sample size for EUR = 100k)")
p
png(file = paste0("./method_compare_result_summary.png"),
    width = 13, height = 8, res = 300,units = "in")
p
dev.off()

#In particular I want the following subcategories
#Single ethnic method: P+T, LDpred
#EUR PRS: EUR P+T PRS, EUR LDpred PRS
#Multi ethnic method: SNPs in EUR P+T + target coefficients, SNPs in EUR P+T + EB coefficients, ME-Bayes, 2DLD, 2DLD-EB

