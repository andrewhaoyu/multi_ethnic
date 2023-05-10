eth = c("AFR","AMR","EAS","EUR","SAS")
herit = c(0.32,0.21,0.16,0.19,0.17)
n.cau = c(14798,9755,7541,8536,7997)
averaged.effect = 0.4/n.cau
data = data.frame(eth,herit,n.cau,averaged.effect)
setwd("/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/LD_simulation_GA")
source("../../code/LD_simulation_large/theme_Publication.R")

p <- ggplot(data,aes(x= eth,y=herit,fill=eth))+
  geom_bar(stat = "identity")+
  theme_Publication()+
  ylab("Heritability")+
  xlab("Ethnic group")+
  #scale_fill_Publication()+
  # scale_color_nejm()+
  ggtitle("Common SNPs heritability")+
  theme(axis.text = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.position = "none")+
  scale_fill_manual(values = c("#386cb0","#EF7E3D","#ffd558","#004953","#7fc97f"))
p 

png(file = paste0("./Heritability_plot.png"),
    width = 7, height = 7, res = 300,units = "in")
print(p)
dev.off()

p <- ggplot(data,aes(x= eth,y=averaged.effect,fill=eth))+
  geom_bar(stat = "identity")+
  theme_Publication()+
  ylab("Average per-SNP heritability")+
  xlab("Ethnic group")+
  # scale_color_nejm()+
  ggtitle("Average per-SNP heritability")+
  theme(axis.text = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.position = "none")+
  scale_fill_manual(values = c("#386cb0","#EF7E3D","#ffd558","#004953","#7fc97f"))
p 

png(file = paste0("./averaged_effect_sizes.png"),
    width = 7, height = 7, res = 300,units = "in")
print(p)
dev.off()
