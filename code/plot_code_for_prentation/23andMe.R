setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/PT_summary/")
eth_group = c("european","african_american",
              "latino","east_asian","south_asian")
eth_name = c("European","African American",
             "Latino","East Asian","South Asian")
eth <- c("EUR","AFR","AMR","EAS","SAS")
bin_trait <- c("any_cvd","depression",
               "iqb.sing_back_musical_note",
               "migraine_diagnosis",
               "morning_person")
con_trait <- c(       "heart_metabolic_disease_burden",
                      "height"
)
trait = c("any_cvd","depression",
          "heart_metabolic_disease_burden",
          "height",
          "iqb.sing_back_musical_note",
          "migraine_diagnosis",
          "morning_person")
trait_name = c("Any CVD","Depression",
               "Heart metabolic disease burden",
               "Height",
               "SBMN",
               "Migraine Diagnosis",
               "Morning Person")

source("/Users/zhangh24/GoogleDrive/multi_ethnic/code/stratch/theme_publication.R")

#load all results
plot.data.list = list()
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(cowplot)
load("/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/prediction_summary.rdata")
prediction.result$method_vec = as.character(prediction.result$method_vec)
prediction.result = prediction.result %>% 
  mutate(method_vec = ifelse(method_vec == "CT-SLEB (two ancestries)",
                             "CT-SLEB",method_vec)) %>% 
  mutate(method_vec = ifelse(method_vec == "Weighted PRS",
                             "Weighted PRS (CT)",method_vec)) 
# prediction.result.table = prediction.result %>% 
#   mutate(sigma2 = ifelse(result<0.5,0,qnorm(result)^2*2)) %>% 
#   mutate(sigma2 = ifelse(trait%in%
#                            c("Heart metabolic disease burden","Height"),NA,sigma2)) %>% 
#   select(eth,trait,Method,result,sigma2)
# write.csv(prediction.result.table,file = "/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/prediction_summary.csv",row.names = F)

uvals = factor(c("CT",
                 "LDpred2",
                 "Best EUR SNP (CT)",
                 "Best EUR SNP (LDpred2)",
                 "Weighted PRS (CT)",
                 "PRS-CSx",
                 "CT-SLEB",
                 "CT-SLEB (five ancestries)"),
               levels= c("CT",
                         "LDpred2",
                         "Best EUR SNP (CT)",
                         "Best EUR SNP (LDpred2)",
                         "Weighted PRS (CT)",
                         "PRS-CSx",
                         "CT-SLEB",
                         "CT-SLEB (five ancestries)"))

n.single = 9
n.EUR = 9

n.multi = 9


single.color =  brewer.pal(n.single, "Blues")[c(4,7)]

EUR.color = brewer.pal(n.EUR, "RdPu")[c(4,7)]

weighted.color = brewer.pal(n.multi, "Greens")[c(3)]

bayes.color = brewer.pal(n.multi, "Oranges")[c(5)]

propose.method = brewer.pal(n.multi, "Purples")[c(4,7)]


colour = c(single.color,EUR.color,weighted.color,
           bayes.color,propose.method)

col_df = tibble(
  colour = c(single.color,EUR.color,weighted.color,
             bayes.color,propose.method),
  method_vec = uvals,
  category = case_when(method_vec%in%c("CT",
                                       "LDpred2") ~ "Single ancestry method",
                       method_vec%in%c("Best EUR SNP (CT)",
                                       "Best EUR SNP (LDpred2)"
                       ) ~ "EUR PRS based method",
                       method_vec%in%c("Weighted PRS (CT)") ~ "Weighted PRS method",
                                       
                       method_vec%in%c("PRS-CSx") ~ "Bayesian method",
                       method_vec%in%c("CT-SLEB",
                         "CT-SLEB (five ancestries)") ~ "Proposed method")) %>%   
  mutate(category = factor(category,levels = c("Single ancestry method",
                                                     "EUR PRS based method",
                                                     "Weighted PRS method",
                                                     "Bayesian method",
                                                     "Proposed method")))
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


legs = lapply(sort(unique(col_df$category)), run_plot)

legs = lapply(legs, getLegend)
p.leg = plot_grid(NULL,legs[[1]],legs[[2]], legs[[3]],legs[[4]],legs[[5]], NULL,align="v",ncol=1,rel_heights=c(1.2,1,1,1,0.8,0.8,0.8,1.2))
print(p.leg)



prediction.result.sub = prediction.result %>% 
  filter(trait %in% c("Any CVD","Depression",
                      "SBMN",
                      "Migraine Diagnosis",
                      "Morning Person")) %>% 
  filter(eth!="European") %>% 
  mutate(eth = factor(eth,
                      levels = c("African American",
                                 "Latino","East Asian","South Asian")))

#create benchmark C+T results for european to plot for other ethnic group
prediction.result.European.CT = prediction.result %>% 
  filter(trait %in% c("Any CVD","Depression",
                      "SBMN",
                      "Migraine Diagnosis",
                      "Morning Person")) %>% 
  filter(eth=="European"&
           method_vec=="CT") %>% 
  select(result,trait)
prediction.result.European.LDpred2 = prediction.result %>% 
  filter(trait %in% c("Any CVD","Depression",
                      "SBMN",
                      "Migraine Diagnosis",
                      "Morning Person")) %>% 
  filter(eth=="European"&
           method_vec=="LDpred2") %>% 
  select(result,trait)

prediction.result.sub = prediction.result.sub %>%
  mutate(sigma2 = ifelse(result<0.5,0,qnorm(result)^2*2))
prediction.result.sub = prediction.result.sub %>% 
  filter(trait=="Any CVD")
prediction.result.European.CT = prediction.result.European.CT %>% 
  filter(trait=="Any CVD")
p.null <- ggplot(prediction.result.sub)+
  geom_bar(aes(x = index,y = result,fill=method_vec),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  ylab("AUC")+
  facet_grid(vars(trait),vars(eth))+
  theme_Publication()+
  coord_cartesian(ylim = c(0.47, 0.67)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values = colour) +
  theme(legend.position = "none")+
  geom_hline(data = prediction.result.European.CT, aes(yintercept =result), linetype = "dashed",color = "red")

print(p.null)
p = plot_grid(p.null,p.leg,nrow=1,rel_widths = c(3,1))
png(filename = "../../presentation_plot/23andme_cvd.png",width=12,height = 8,units = "in",res = 300)
print(p)
dev.off()



prediction.result.European = prediction.result %>% 
  filter(trait %in% c("Heart metabolic disease burden")==T&
           method_vec=="CT") %>% 
  filter(eth=="European") %>% 
  select(result,trait)
prediction.result.sub = prediction.result %>% 
  filter(trait %in% c("Heart metabolic disease burden")) %>% 
  filter(eth!="European") %>% 
  mutate(eth = factor(eth,
                      levels = c("African American",
                                 "Latino","East Asian","South Asian")))

p.null <- ggplot(prediction.result.sub)+
  geom_bar(aes(x = index,y = result,fill=method_vec),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  ylab("R2")+
  facet_grid(vars(trait),vars(eth),scales = "free")+
  theme_Publication()+
  #coord_cartesian(ylim = c(0.47, 0.65)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values = colour) +
  theme(legend.position = "none")+
  geom_hline(data = prediction.result.European, aes(yintercept = result), linetype = "dashed",color = "red")



p = plot_grid(p.null,p.leg,nrow=1,rel_widths = c(3,1))
#theme(strip.text.y = element_blank())
#library(RColorBrewer)
png(filename = "../../presentation_plot/23andme_heart.png",width=12,height = 8,units = "in",res = 300)
print(p)
dev.off()




prediction.result.European = prediction.result %>% 
  filter(trait %in% c("Height")==T&
           method_vec=="CT") %>% 
  filter(eth=="European") %>% 
  select(result,trait)
prediction.result.sub = prediction.result %>% 
  filter(trait %in% c("Height")) %>% 
  filter(eth!="European") %>% 
  mutate(eth = factor(eth,
                      levels = c("African American",
                                 "Latino","East Asian","South Asian")))

p.null <- ggplot(prediction.result.sub)+
  geom_bar(aes(x = index,y = result,fill=method_vec),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  ylab("R2")+
  facet_grid(vars(trait),vars(eth),scales = "free")+
  theme_Publication()+
  #coord_cartesian(ylim = c(0.47, 0.65)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values = colour) +
  theme(legend.position = "none")+
  geom_hline(data = prediction.result.European, aes(yintercept = result), linetype = "dashed",color = "red")



p = plot_grid(p.null,p.leg,nrow=1,rel_widths = c(3,1))
#theme(strip.text.y = element_blank())
#library(RColorBrewer)
png(filename = "../../presentation_plot/23andme_height.png",width=12,height = 8,units = "in",res = 300)
print(p)
dev.off()

