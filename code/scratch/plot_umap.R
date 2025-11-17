
# input files: Y.full_score.gz, dataUMAP.txt
setwd("./data/UMAP_plot/")
dataScore <- read.table("Y.full_score.gz", header = TRUE)
#dim(dataScore)
#head(dataScore[, 1:7])
#hist(dataScore$norm_score, 100)

dataUMAP <- read.table("dataUMAP.txt", header = TRUE, row.names = 1)
#str(dataUMAP)
#plot(dataUMAP)

colorMap <- colorRampPalette(c("blue", "white", "red"))(100)
tempColor <- colorMap[round((dataScore$norm_score + 5) * 10)]

plot(dataUMAP, col = tempColor, pch = 20)

library(RColorBrewer)
library(viridis)
library(scico)
data = cbind(dataUMAP,dataScore)
#mid <- mean(data$norm_score)
library(ggplot2)
p = ggplot(data = data) + 
  geom_point(aes(x = UMAP1, y = UMAP2, col = norm_score)) +
  scale_color_gradientn(colors = rev(brewer.pal(10, "RdBu")))+
  theme_Publication()+
  labs(color = "Disease Score")
