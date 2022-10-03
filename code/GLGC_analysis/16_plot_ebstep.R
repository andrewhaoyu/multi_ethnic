setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/GLGC/analysis_result/")
source("/Users/zhangh24/GoogleDrive/multi_ethnic//code/LD_simulation_large/theme_Publication.R")
library(ggplot2)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(data.table)
#library(RColorBrewer)
load("ct_sleb_ebfirst.rdata")
ct_sleb_test= final_result %>% 
  mutate(method = "Method 1")
load("ct_sleb_all.rdata")
ct_sleb_all = final_result

prediction.result <- rbind(ct_sleb_all,
                           ct_sleb_test) %>% 
  rename(result = r2)

prediction.result$index = rep("1",nrow(prediction.result))
prediction.result = prediction.result %>% 
  mutate(trait = 
           case_when(trait == "HDL" ~ "HDL",
                     trait == "LDL" ~ "LDL",
                     trait == "logTG" ~ "logTG",
                     trait == "TC" ~ "TC"))

prediction.result.sub = prediction.result %>% 
  filter(eth%in%c("EUR","AMR") == F) 



p.null <- ggplot(prediction.result.sub)+
  geom_bar(aes(x = index,y = result,fill=method),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  ylab(expression(bold(R^2)))+
  facet_grid(vars(trait),vars(eth),scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  #coord_cartesian(ylim = c(0.47, 0.67)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
print(p.null)
#scale_fill_manual(values = colour) 




n.single = 9
n.EUR = 9

n.multi = 9


single.color =  brewer.pal(n.single, "Blues")[c(4,7)]

EUR.color = brewer.pal(n.EUR, "RdPu")[c(4,7)]

weighted.color = brewer.pal(n.multi, "Greens")[c(3,5,7)]

bayes.color = brewer.pal(n.multi, "Oranges")[c(4,7)]

propose.method = brewer.pal(n.multi, "Purples")[c(3,5,7)]



colour = c(single.color,EUR.color,weighted.color,
           bayes.color,propose.method)

col_df = tibble(
  colour = c(single.color,EUR.color,weighted.color,
             bayes.color,propose.method),
  method_vec = uvals,
  category = case_when(method_vec%in%c("CT",
                                       "LDpred2") ~ "Single ancestry method",
                       method_vec%in%c("Best EUR SNP (CT)",
                                       "Best EUR PRS (LDpred2)"
                       ) ~ "EUR PRS based method",
                       method_vec%in%c("Weighted PRS (CT)",
                                       "Weighted PRS (LDpred2)",
                                       "PolyPred+"
                       ) ~ "Weighted PRS method",
                       method_vec%in%c("XPASS",
                                       "PRS-CSx") ~"Bayesian method",
                       method_vec%in%c("CT-SLEB", "MEBayes", "MELasso") ~ "Proposed method")
) %>% 
  mutate(category = factor(category,levels = c("Single ancestry method",
                                               "EUR PRS based method",
                                               "Weighted PRS method",
                                               "Bayesian method",
                                               "Proposed method")))

prediction.result.European =prediction.result %>% 
  filter(eth=="EUR"&
           method_vec%in%c("CT","LDpred2")) %>% 
  group_by(trait) %>% 
  summarise(newresult = max(result))

prediction.result = prediction.result %>% 
  left_join(col_df)
getLegend <- function(p) {
  g <- ggplotGrob(p)
  k <- which(g$layout$name=="guide-box")
  g$grobs[[k]]
}

run_plot = function(filler, values) {
  values = col_df %>% 
    filter(category %in% filler)
  labels = values %>% 
    pull(method_vec)
  values = values %>% pull(colour)
  names(values) = labels
  ggplot(
    prediction.result %>% 
      filter(category %in% filler),
    aes(x= index,y=result,
        group=method_vec))+
    geom_bar(aes(fill = method_vec),
             stat="identity",
             position = position_dodge())+
    #geom_point(aes(color=method_vec))+
    theme_Publication()+
    ylab("R2")+
    xlab("Sample Size")+
    labs(fill = filler)+
    scale_fill_manual(values = values)
}

library(cowplot)
legs = lapply(sort(unique(col_df$category)), run_plot)

legs = lapply(legs, getLegend)
p.leg = plot_grid(NULL,legs[[1]],legs[[2]], legs[[3]],legs[[4]],legs[[5]], NULL,align="v",ncol=1,rel_heights=c(1.2,1,1,1,1,1,1,1.2))
print(p.leg)





p.null <- ggplot(prediction.result.sub)+
  geom_bar(aes(x = index,y = result,fill=method_vec),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  ylab(expression(bold(R^2)))+
  facet_grid(vars(trait),vars(eth),scales = "free")+
  theme_Publication()+
  #coord_cartesian(ylim = c(0.47, 0.67)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values = colour) +
  theme(legend.position = "none") +
  geom_hline(data = prediction.result.European, aes(yintercept = newresult), linetype = "dashed",color = "red")

print(p.null)
p = plot_grid(p.null,p.leg,nrow=1,rel_widths = c(3,1))
png(filename = "./GLGC_result_all.png",width=17,height = 12,units = "in",res = 300)
print(p)
dev.off()
