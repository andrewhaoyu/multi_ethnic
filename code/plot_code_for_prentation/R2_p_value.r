setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01,0.5)
R2 = c(0.04,0.042,0.045,0.046,0.052,0.056,0.054,0.050,0.048,0.042,0.039,0.032)
plot.data = data.frame(x = log10(pthres),R2) 
library(tidyverse)
library(ggplot2)

p = ggplot(plot.data,aes(x,R2))+
  geom_point(size = 2)+
  geom_line(size = 1.5)+theme_Publication()+
  #theme(axis.text.y=element_blank())+
  ylab("R2")+
  xlab("log10(P)")+
  ggtitle(NULL)+
  theme(axis.text = element_text(size = rel(1.8)),
        axis.title = element_text(face = "bold",size = rel(1.8)))
png(file = paste0("./result/presentation_plot/R2_pvalue.png"),width = 8, height = 6, units = "in",res = 300)
print(p)
dev.off()







y = dgamma(x,shape = 6, rate = 1.5)
plot(x,y)
risk = seq(0,60,by = 60/(length(x)-1))
library(dplyr)
plot.data = data.frame(risk,y) %>% 
  mutate(risk_category = case_when(risk>25~"High risk",
                                   risk>15~"Moderate risk",
                                   TRUE~"Low risk")) %>% 
  mutate(risk_category = factor(risk_category,
                                levels = c("Low risk",
                                           "Moderate risk",
                                           "High risk")))


p = ggplot(plot.data,aes(risk,y))+
  geom_line()+theme_Publication()+
  theme(axis.text.y=element_blank())+
  ylab("Frequncy")+
  xlab("Lifetime absolute risk (%)")+
  geom_area(aes(fill =risk_category ))+
  scale_fill_brewer(palette = "OrRd")+
  geom_vline(xintercept = 15,linetype = "dashed")+
  geom_vline(xintercept = 25,linetype = "dashed")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle("Absolute risk")
png(file = paste0("./result/presentation_plot/AFR_risk_stratition.png"),width = 12, height = 6, units = "in",res = 300)
print(p)
dev.off()

