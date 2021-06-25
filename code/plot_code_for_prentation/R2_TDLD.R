#demontration plot
#load AMR TDLD results for presentation
#l = 2
#m = 1
#i1 = 1

setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result")
load("AMR_TDLD_for_plot.rdata")



source("../code/independent_SNPs_simulation/regression_prediciton/plot_theme.R")
library(tidyverse)


plot.data = result.data %>%
  filter(r2_ind_vec == 3&
           wc_ind_vec==1) %>%
  rename(r2 = r2.vec.test) %>% 
  select(r2,pthres_vec1,pthres_vec2)

#add p-value as 1.0 for plot
pthres_vec1 = c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)
pthres_vec2 = c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01)

r2 = runif(length(pthres_vec1),0.9*0.01254517,0.9*0.01400249)

plot.data.sub = data.frame(pthres_vec1,pthres_vec2,r2)

plot.data = rbind(plot.data,plot.data.sub)

plot.data = plot.data %>%
  mutate(peur = as.factor(-round(log10(pthres_vec1),1)),
         ptar = as.factor(-round(log10(pthres_vec2),1))) 

p= ggplot(data = plot.data, aes(x=peur, y=ptar)) +
  geom_tile(aes(fill=r2),colour = "white")+
  scale_fill_gradient2(low = "grey94",
                       high = "dodgerblue4") +
  fte_theme()+
  labs(x = paste0("EUR -log10(P-value)"),
       y = paste0("Tar -log10(P-value)")) +
  ggtitle(NULL)+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(fill = "R2")+
  theme(axis.text.x=element_text(size=12,color="black",face="bold")) +
  theme(axis.text.y=element_text(size=12,color="black",face="bold")) +
  theme(legend.text = element_text(size=12,color="black")) +
  theme(legend.title = element_text(size=12,color="black",face="bold"))

png(filename="./presentation_plot/TDLD_demonstration.png",
    width = 8, height = 6, res =600,units = "in")
print(p)
dev.off()




