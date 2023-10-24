#plot PCA results in hapgen2
library(data.table)
library(ggplot2)
PCs = fread("./data/pca_005.eigenvec")
source("./code/stratch/theme_publication.R")
PCs$ETH<-substr(PCs$IID,1,3)
colour = c("#1f77b4",
           "#9467bd",
           "#2ca02c",
  "#ff7f0e",
           "#d62728")
p1 = ggplot(PCs,aes(x=PC1,y=PC2,color=as.factor(ETH)))+
  geom_point()+
  theme_Publication()+
  labs(color = "Ancestries")+
  scale_fill_manual(values = colour) 
png(file = paste0("/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/LD_simulation_GA/LD_stack/simulated_data_PC_plot.png"),
    width = 10, height = 8, res = 300,units = "in")
print(p1)
dev.off()
