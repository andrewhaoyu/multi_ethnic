library(reshape2)
library(ggplot2)
library(RColorBrewer)
A <- matrix(rbinom(8*10,2,0.33), 
            byrow = TRUE, nrow = 8, ncol = 10)
colnames(A) <- letters[1:10]
rownames(A) <- LETTERS[1:8]
longData <- melt(A)
fill_color = brewer.pal(9, "Greens")[c(1,5,8)]


p = ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = as.character(value)), color = "black",size = 0.7) +  # Add black borders
  scale_fill_manual(values = fill_color)+
  labs(x = NULL, y = NULL, title = "Matrix") +  # Remove axis labels
  theme_Publication() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",  # Remove legend
    plot.title = element_blank(),
    axis.line = element_blank()  # Remove axis lines
  )
png(file = "./generate_plot/result/")
