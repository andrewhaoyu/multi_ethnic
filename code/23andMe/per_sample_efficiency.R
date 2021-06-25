setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/PT_summary/summary")
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
method_vec = c("PT")
method_name = c("P + T")
source("/Users/zhangh24/GoogleDrive/multi_ethnic/code/stratch/theme_publication.R")
library(tidyverse)
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
method = "PT"
#load all results
plot.data.list = list()
sample_size = read.csv("/Users/zhangh24/GoogleDrive/multi_ethnic/data/23_sample_size.csv")

library(data.table)
temp = 1
for(i1 in 1:length(method_vec)){
  for(l in 1:7){
    for(i in 1:5){
      result = read.csv(paste0(eth_group[i],"_",trait[l],"_",method_vec[i1]))
      
      plot.data = data.frame(result = as.numeric(max(result)),
                             eth = eth_name[i],
                             trait = trait_name[l],
                             Method = method_name[i1])
      
      plot.data.list[[temp]] = plot.data
      temp = temp+1
    }
  }
}
plot.data = rbindlist(plot.data.list)

plot.data$index = rep("1",nrow(plot.data))

plot.data$eth = factor(plot.data$eth,
                       levels = c("European","African American",
                                  "Latino","East Asian","South Asian"))
plot.data$N = sample_size$N_effect

plot.data = plot.data %>% 
  mutate(efficiency = case_when(trait %in%bin_trait==T ~
                                  2*qnorm(result)^2/N,
         TRUE ~ result/N))

plot.data.sub = plot.data %>% 
  filter(trait %in% c("Any CVD","Depression",
                      "SBMN",
                      "Migraine Diagnosis",
                      "Morning Person")) 







p = ggplot(plot.data.sub)+
  geom_bar(aes(x = index,y =  efficiency,fill=eth),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  ylab("Efficiency")+
  facet_wrap(vars(trait))+
  theme_Publication()+
  scale_fill_Publication()+
 # coord_cartesian(ylim = c(0.5, 0.64)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#library(RColorBrewer)

png(filename = "../../comparasion_plot/efficienty_bin.png",width=12,height = 8,units = "in",res = 300)
print(p)
dev.off()

plot.data.sub = plot.data %>% 
  filter(trait %in% c("Heart metabolic disease burden","Height")) 


p = ggplot(plot.data.sub)+
  geom_bar(aes(x = index,y = efficiency,fill=eth),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  ylab("Efficiency")+
  facet_wrap(vars(trait))+
  #facet_wrap(vars(eth))+
  theme_Publication()+
  scale_fill_Publication()+
  #coord_cartesian(ylim = c(0.5, 0.64)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#theme(strip.text.y = element_blank())
#library(RColorBrewer)
png(filename = "../../comparasion_plot/con_comparasion.png",width=10,height = 8,units = "in",res = 300)
print(p)
dev.off()





sample_size = read.csv("/Users/zhangh24/GoogleDrive/multi_ethnic/data/23_sample_size.csv",stringsAsFactors = F)
sample_size = sample_size %>% 
  mutate(N_case = as.numeric(N_case)) %>% 
  mutate(N_case = ifelse(is.na(N_case),0,N_case),
         total = N_control + N_case)
sample_size %>% group_by(eth) %>% 
  summarise(mean(total)/0.70)
