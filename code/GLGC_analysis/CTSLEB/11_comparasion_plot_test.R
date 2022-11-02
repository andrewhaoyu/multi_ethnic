#goal: plot the five ancetries eb results with two ans posterior and five ans posterior
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
load("ct_sleb_all_test.rdata")
# final_result %>% filter(eth == "SAS"& trait == "LDL")
ct_sleb_test1= final_result 

load("ct_sleb_all.rdata")
ct_sleb_all = final_result

prediction.result <- rbind(ct_sleb_all,
                           ct_sleb_test1) %>% 
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
  scale_fill_Publication()+
  #coord_cartesian(ylim = c(0.47, 0.67)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#print(p.null)
png(filename = "./GLGC_large_effect_snp_test.png",width=17,height = 12,units = "in",res = 300)
print(p.null)
dev.off()
