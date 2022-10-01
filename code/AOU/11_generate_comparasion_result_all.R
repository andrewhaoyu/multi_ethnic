setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/AOU/analysis_result/")
source("/Users/zhangh24/GoogleDrive/multi_ethnic//code/LD_simulation_large/theme_Publication.R")
library(ggplot2)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(data.table)
#library(RColorBrewer)
load("PT.rdata")
LD.clump.result = final_result
load("best_eur.rdata")
eursnp.result = final_result
eursnp.result = eursnp.result %>% 
  mutate(method = "Best EUR SNP (CT)")
load("weighted_prs.rdata")
weightedprs.result = final_result %>% 
  mutate(method = "Weighted PRS (CT)")
load("polypred.rdata")
polypred.result = final_result %>% 
  mutate(method = "PolyPred+")
load("xpass.rdata")
xpass.result = final_result


load("prscsx_all.rdata")
prscsx.all.result = final_result %>% mutate(method = "PRS-CSx")

load("ct_sleb_all.rdata")
ct.sleb.all = final_result %>% mutate(method = "CT-SLEB")
load("allofus.jin.Rdata")
ldpred2.result = res %>% 
  filter(Method%in% c("LDpred2", "EUR LDpred2", "Weighted LDpred2","MEBayes SL")) %>% 
  mutate(method = case_when(
    Method == "EUR LDpred2" ~ "Best EUR PRS (LDpred2)",
    Method == "Weighted LDpred2" ~ "Weighted PRS (LDpred2)",
    Method == "LDpred2" ~ "LDpred2",
    Method == "MEBayes SL" ~ "MEBayes"
  )) %>% 
  rename(trait = Trait,
         eth = Race,
         r2 = R2) %>% 
  filter(trait!= "nonHDL") %>% 
  mutate(trait = as.character(trait)) %>% 
  select(eth, trait, method, r2) 

lasso = fread("allofus_jingning.txt")
lasso.result = lasso %>% 
  filter(method_vec%in% c("multi-ethnic lasso (super learning)")) %>% 
  mutate(method = case_when(
    method_vec == "multi-ethnic lasso (super learning)" ~ "MELasso"
  )) %>% 
  rename(trait = t_vec,
         eth = eth.vec,
         r2 = r2.vec) %>% 
  filter(trait!= "nonHDL"  ) %>% 
  select(eth, trait, method, r2) 





prediction.result <- rbind(LD.clump.result,
                           eursnp.result,
                           ldpred2.result,
                           weightedprs.result,
                           polypred.result,
                           xpass.result,
                           prscsx.all.result,
                           ct.sleb.all,
                           lasso.result) %>% 
  rename(result = r2)

prediction.result = prediction.result %>% 
  mutate(method_vec = method) %>% 
  mutate(
    method_vec = factor(method_vec,
                        levels = c("CT",
                                   "LDpred2",
                                   "Best EUR SNP (CT)",
                                   "Best EUR PRS (LDpred2)",
                                   "Weighted PRS (CT)",
                                   "Weighted PRS (LDpred2)",
                                   "PolyPred+", 
                                   "XPASS", 
                                   "PRS-CSx",
                                   "CT-SLEB",
                                   "MEBayes",
                                   "MELasso"
                        )))


uvals = factor(c("CT",
                 "LDpred2",
                 "Best EUR SNP (CT)",
                 "Best EUR PRS (LDpred2)",
                 "Weighted PRS (CT)",
                 "Weighted PRS (LDpred2)",
                 "PolyPred+", 
                 "XPASS", 
                 "PRS-CSx",
                 "CT-SLEB",
                 "MEBayes",
                 "MELasso"
)
,levels = c("CT",
            "LDpred2",
            "Best EUR SNP (CT)",
            "Best EUR PRS (LDpred2)",
            "Weighted PRS (CT)",
            "Weighted PRS (LDpred2)",
            "PolyPred+", 
            "XPASS", 
            "PRS-CSx",
            "CT-SLEB",
            "MEBayes",
            "MELasso"
))

prediction.result$index = rep("1",nrow(prediction.result))
prediction.result = prediction.result %>% 
  mutate(trait = 
           case_when(trait == "height" ~ "Height",
                     trait == "bmi" ~ "BMI"))
save(prediction.result,file = "aou.prediction.result.summary.all.rdata")
prediction.result.sub = prediction.result %>% 
  filter(eth%in%c("EUR","AMR") == F) 
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
png(filename = "./aou_result_all.png",width=17,height = 12,units = "in",res = 300)
print(p)
dev.off()
