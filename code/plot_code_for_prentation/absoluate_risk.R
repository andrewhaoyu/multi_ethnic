x = seq(0,20,0.005)
y = dgamma(x,shape = 2, rate = 0.5)
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
library(ggplot2)

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
  ggtitle("Absolute risk")+
   guides(fill=guide_legend(title="Risk category"))
png(file = paste0("./result/presentation_plot/EUR_risk_stratition.png"),width = 12, height = 6, units = "in",res = 300)
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

