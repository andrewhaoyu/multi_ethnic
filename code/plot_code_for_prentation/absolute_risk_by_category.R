risk = seq(0,20,0.005)
high_mean = 2.1
high_sd = 0.9
low_mean = 1.5
low_sd = 0.25
median_mean = 1.7
median_sd = DifferenceCal(low_mean,high_sd,median_mean)
LogNormalMeanEst = function(logmean, logsd){
  exp(logmean + logsd^2/2)
}

DifferenceCal = function(low_mean,high_sd,median_mean){
  sqrt(2*(log(LogNormalMeanEst(low_mean,high_sd))-median_mean))
}

LogNormalMeanEst(high_mean,low_sd)
LogNormalMeanEst(low_mean,high_sd)
LogNormalMeanEst(median_mean,median_sd)
P1_mean = LogNormalMeanEst(high_mean,low_sd)
P2_mean = LogNormalMeanEst(low_mean,high_sd)
P3_mean = LogNormalMeanEst(median_mean,median_sd)
y = dlnorm(risk,mean = high_mean, sd = low_sd)
mean_est_high = mean(rlnorm(20000,meanlog = high_mean, sdlog = high_sd))
print(mean_est_high)

library(dplyr)
library(RColorBrewer)
plot.data = data.frame(risk,y) %>% 
  mutate(risk_category = case_when(risk>10~"High risk",
                                   # risk>15~"Moderate risk",
                                   TRUE~"Low risk")) %>% 
  mutate(risk_category = factor(risk_category,
                                levels = c("Low risk",
                                           #            "Moderate risk",
                                           "High risk")))

library(ggplot2)
colour = brewer.pal(9, "OrRd")[c(3,7)]
p1 = ggplot(plot.data,aes(risk,y))+
  theme_Publication()+
  theme(axis.text.y=element_blank())+
  ylab(NULL)+
  xlab(NULL)+
  geom_area(aes(fill =risk_category ))+
  scale_fill_manual(values = colour) +
  # geom_vline(xintercept = 15,linetype = "dashed")+
  # geom_vline(xintercept = 25,linetype = "dashed")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  
  ggtitle(NULL)+
  guides(fill=guide_legend(title="Risk category"))+
  geom_line(size = 1)+
  geom_point(aes(x = P1_mean, y = 0), colour = "red", shape = 8, size = 4)
#   theme(
#           axis.line = element_line(colour="black",size=1)
#   )
  
  
#print(p1)

risk = seq(0,20,0.005)
y = dlnorm(risk,mean =low_mean, sd = high_sd)
#y = dgamma(x,shape =2, rate = 0.5)

library(dplyr)
library(RColorBrewer)
plot.data = data.frame(risk,y) %>% 
  mutate(risk_category = case_when(risk>10~"High risk",
                                   # risk>15~"Moderate risk",
                                   TRUE~"Low risk")) %>% 
  mutate(risk_category = factor(risk_category,
                                levels = c("Low risk",
                                           #            "Moderate risk",
                                           "High risk")))
library(ggplot2)
colour = brewer.pal(9, "OrRd")[c(3,7)]
p2 = ggplot(plot.data,aes(risk,y))+
  theme_Publication()+
  theme(axis.text.y=element_blank())+
  ylab(NULL)+
  xlab(NULL)+
  geom_area(aes(fill =risk_category ))+
  scale_fill_manual(values = colour) +
  # geom_vline(xintercept = 15,linetype = "dashed")+
  # geom_vline(xintercept = 25,linetype = "dashed")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  ggtitle(NULL)+
  guides(fill=guide_legend(title="Risk category"))+
  geom_line(size = 1)+
  geom_point(aes(x = P2_mean, y = 0), colour = "red", shape = 8, size = 4)
#print(p2)



risk = seq(0,20,0.005)
y = dlnorm(risk,mean =median_mean, sd = median_sd)
#y = dgamma(x,shape =2, rate = 0.5)

library(dplyr)
library(RColorBrewer)
plot.data = data.frame(risk,y) %>% 
  mutate(risk_category = case_when(risk>10~"High risk",
                                   # risk>15~"Moderate risk",
                                   TRUE~"Low risk")) %>% 
  mutate(risk_category = factor(risk_category,
                                levels = c("Low risk",
                                           #            "Moderate risk",
                                           "High risk")))
library(ggplot2)
colour = brewer.pal(9, "OrRd")[c(3,7)]
p3 = ggplot(plot.data,aes(risk,y))+
  theme_Publication()+
  theme(axis.text.y=element_blank())+
  ylab(NULL)+
  xlab("Absolute risk")+
  geom_area(aes(fill =risk_category ))+
  scale_fill_manual(values = colour) +
  # geom_vline(xintercept = 15,linetype = "dashed")+
  # geom_vline(xintercept = 25,linetype = "dashed")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  geom_line(size = 1)+
  ggtitle(NULL)+
  guides(fill=guide_legend(title="Risk category"))+
  geom_point(aes(x = P2_mean, y = 0), colour = "red", shape = 8, size = 4)
#print(p3)

library(cowplot)
p = plot_grid(p1,p2,p3,ncol = 1)
print(p)
png(file = paste0("./result/presentation_plot/absolute_demonstration_figure.png"),width = 14, height = 16, units = "in",res = 300)
print(p)
dev.off()
