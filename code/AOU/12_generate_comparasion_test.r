setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/AOU/analysis_result/")
source("/Users/zhangh24/GoogleDrive/multi_ethnic//code/LD_simulation_large/theme_Publication.R")
library(ggplot2)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(data.table)

load("ct_sleb_all.rdata")
ct.sleb.all = final_result %>% mutate(method = "CT-SLEB")
load("ct_sleb_five_cut_test.rdata")
ct.sleb.all.test = final_result %>% 
  filter(method == "CT-SLEB (pthres =1e-10)") %>% 
  mutate(method = "CT-SLEB (large effect SNPs not in EB)")

prediction.result <- rbind(ct.sleb.all,
                           ct.sleb.all.test) %>% 
  rename(result = r2)

prediction.result$index = rep("1",nrow(prediction.result))
prediction.result = prediction.result %>% 
  mutate(trait = 
           case_when(trait == "height" ~ "Height",
                     trait == "bmi" ~ "BMI"))

prediction.result.sub = prediction.result %>% 
  filter(eth%in%c("EUR","AMR") == F) 

p.null <- ggplot(prediction.result.sub)+
  geom_bar(aes(x = index,y = result,fill=method),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  ylab(expression(bold(R^2)))+
  facet_grid(vars(trait),vars(eth),scales = "free")+
  scale_fill_Publication()+
  #coord_cartesian(ylim = c(0.47, 0.67)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
png(filename = "./AOU_large_effect_snp_test.png",width=17,height = 12,units = "in",res = 300)
print(p.null)
dev.off()

  