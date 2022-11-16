setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA/LD_stack/")
source("../../../code/LD_simulation_large/theme_Publication.R")

label <- paste0("X", 1:6)
mean  <- c(1.29,0.76,2.43,1.68,1.22,1.7) 
lower <- c(1.29,0.76,2.43,1.68,1.22,1.7)
upper <- c(1.29,0.76,2.43,1.68,1.22,1.7)

df <- data.frame(label, mean, lower, upper)

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label))

library(ggplot2)
fp <- ggplot(data=df, aes(x=label, y=mean)) +
  geom_point( color = "red", size = 5) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Label") + ylab("Mean (95% CI)") +
  theme_Publication()+ # use a white background  
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),panel.border = element_rect(fill = NA, 
                                                                colour = "black",size=3),
        axis.line = element_line(colour="black",size=0.7)) 
print(fp)

data_table <- ggplot(data = df, aes(y = label)) +
  #geom_hline(aes(yintercept = label), size = 7) +
  geom_text(aes(x = 0, label = label), hjust = 0) +
  #geom_text(aes(x = 5, label = eventnum)) +
  #geom_text(aes(x = 7, label = arr), hjust = 1) +
  scale_colour_identity() +
  theme_void() + 
  theme(plot.margin = margin(5, 0, 35, 0))
grid.arrange(data_table,fp, ncol = 2)
