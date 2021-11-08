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
method_vec = c("PT","BESTEUR","weighted_PRS","tdld","TDLD_EB","TDLD_SLEB","TDLD_SLEBall")
method_name = c("C+T","Best EUR SNP (C+T)",
                "Weighted PRS",
                "TDLD",
                "TDLD-EB",
                "TDLD-SLEB",
                "TDLD-SLEB (all ethnics)")
besteur_methodname = c("Best EUR SNP (C+T)",
                       "Best EUR SNP + target coefficients (C+T)",
                       "Best EUR SNP + EB coefficients (C+T)")
source("/Users/zhangh24/GoogleDrive/multi_ethnic/code/stratch/theme_publication.R")
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
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
temp = 1

for(l in 1:7){
  for(i in 2:5){
    for(i1 in 1:length(method_vec)){
      if(method_vec[i1]=="BESTEUR"){
        #BEST EUR doens't need tuning parameters
        result = read.csv(paste0("./testing_summary/",eth_group[i],"_",trait[l],"_",method_vec[i1]))
        #three different BEST EUR methods
        for(k in 1:3){
          plot.data = data.frame(result = result[k],
                                 eth = eth_name[i],
                                 trait = trait_name[l],
                                 method_vec = besteur_methodname[k])
          colnames(plot.data)[1] = "result"
          plot.data.list[[temp]] = plot.data
          
          temp = temp+1
          
        }
        
      }else if(method_vec[i1]%in%c("weighted_PRS","TDLD_SLEB","TDLD_SLEBall")){
        #weighted PRS, TDLD-SLEB and TDLD-SLEB all ethnic only need validation
        result = read.csv(paste0("./summary_",method_vec[i1],"/",eth_group[i],"_",trait[l]))
        plot.data = data.frame(result = result,
                               eth = eth_name[i],
                               trait = trait_name[l],
                               method_vec = method_name[i1])
        colnames(plot.data)[1] = "result"
        plot.data.list[[temp]] = plot.data
        temp = temp+1
      }else{
        #the other methods need to separate tuning and validation dataset
        result = read.csv(paste0("./validation_summary/",eth_group[i],"_",trait[l],"_",method_vec[i1]))
        idx.max = which.max(result)
        result = read.csv(paste0("./testing_summary/",eth_group[i],"_",trait[l],"_",method_vec[i1]))
        #result2 = read.csv(paste0("./validation_summary/",eth_group[i],"_",trait[l],"_",method_vec[i1]))
        plot.data = data.frame(result = result[idx.max],
                               eth = eth_name[i],
                               trait = trait_name[l],
                               method_vec = method_name[i1])
        colnames(plot.data)[1] = "result"
        plot.data.list[[temp]] = plot.data
        temp = temp+1
      }
      
      
    }
  }
}
prediction.result = rbindlist(plot.data.list)
prediction.result$method_vec = factor(prediction.result$method_vec,
                                      levels = c("C+T","Best EUR SNP (C+T)",
                                                 "Best EUR SNP + target coefficients (C+T)",
                                                 "Best EUR SNP + EB coefficients (C+T)",
                                                 "Weighted PRS",
                                                 "TDLD",
                                                 "TDLDEB",
                                                 "TDLD-EB",
                                                 "TDLD-SLEB",
                                                 "TDLD-SLEB (all ethnics)"))
prediction.result$Method = prediction.result$method_vec
prediction.result$index = rep("1",nrow(prediction.result))
prediction.result$eth = factor(prediction.result$eth,
                               levels = c("European","African American",
                                          "Latino","East Asian","South Asian"))
uvals = unique(prediction.result$method_vec)

n.single = 9


single.color =  brewer.pal(n.single, "Blues")[c(5)]
n.EUR = 9


EUR.color = brewer.pal(n.EUR, "Greens")[c(4,5,6)]


n.multi = 9
multi.color = brewer.pal(n.multi, "Oranges")[c(3:7)]
colour = c(single.color,EUR.color,multi.color)
col_df = tibble(
  colour = c(single.color,EUR.color,multi.color),
  method_vec = uvals,
  category = case_when(method_vec%in%c("C+T") ~ "Single ethnic method",
                       method_vec%in%c("Best EUR SNP (C+T)",
                                       "Best EUR SNP + target coefficients (C+T)",
                                       "Best EUR SNP + EB coefficients (C+T)"
                       ) ~ "EUR PRS based method",
                       method_vec%in%c("Weighted PRS",
                                       "TDLD",
                                       "TDLD-EB",
                                       "TDLD-SLEB",
                                       "TDLD-SLEB (all ethnics)") ~ "Multi ethnic method")
) %>%   mutate(category = factor(category,levels = c("Single ethnic method",
                                                     "EUR PRS based method",
                                                     "Multi ethnic method")))

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
    geom_bar(aes(fill = Method),
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
p.leg = plot_grid(NULL,NULL,legs[[1]],legs[[2]], legs[[3]],NULL,align="v",ncol=1,rel_heights=c(1,1,0.7,1,1.2,1))
print(p.leg)



prediction.result.sub = prediction.result %>% 
  filter(trait %in% c("Any CVD","Depression",
                      "SBMN",
                      "Migraine Diagnosis",
                      "Morning Person")) 

p.null <- ggplot(prediction.result.sub)+
  geom_bar(aes(x = index,y = result,fill=method_vec),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  ylab("AUC")+
  facet_grid(vars(trait),vars(eth))+
  theme_Publication()+
  coord_cartesian(ylim = c(0.47, 0.65)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values = colour) +
  theme(legend.position = "none")


p = plot_grid(p.null,p.leg,nrow=1,rel_widths = c(3,1))
png(filename = "../comparasion_plot/bin_comparasion.png",width=12,height = 8,units = "in",res = 300)
print(p)
dev.off()

#convert AUC to liability scale variance
#AUC = pnorm(sqrt(sigma^2/2))
#sigma^2 = qnorm(AUC)^2*2
sigma2toauc = function(x){
  ifelse(x==0,0.50,round(pnorm(0.5*sqrt(x)),2))
}
prediction.result.sub = prediction.result.sub %>%
  mutate(sigma2 = ifelse(result<0.5,0,qnorm(result)^2*2))

p.null <- ggplot(prediction.result.sub)+
  geom_bar(aes(x = index,y = sigma2,fill=method_vec),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  scale_y_continuous(labels = sigma2toauc)+
  ylab("AUC")+
  facet_grid(vars(trait),vars(eth))+
  theme_Publication()+
  #coord_cartesian(ylim = c(0.47, 0.65)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values = colour) +
  theme(legend.position = "none")

p = plot_grid(p.null,p.leg,nrow=1,rel_widths = c(3,1))
png(filename = "../comparasion_plot/bin_comparasion_r2.png",width=12,height = 8,units = "in",res = 300)
print(p)
dev.off()




prediction.result.sub = prediction.result %>% 
  filter(trait %in% c("Heart metabolic disease burden","Height")) 


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
  theme(legend.position = "none")


p = plot_grid(p.null,p.leg,nrow=1,rel_widths = c(3,1))
#theme(strip.text.y = element_blank())
#library(RColorBrewer)
png(filename = "../comparasion_plot/con_comparasion.png",width=12,height = 8,units = "in",res = 300)
print(p)
dev.off()





#change AUC to R2
prediction.result.sub = prediction.result %>% 
  mutate(R2 = ifelse(trait %in% c("Heart metabolic disease burden","Height"),
                     result,
                     ifelse(result<0.5,0,qnorm(result)^2*2))) %>% 
  filter(method_vec%in%c("Weighted PRS","TDLD-SLEB"))

#generate relative R2 improvment
trait = rep("c",1000)
eth = rep("c",1000)
relative_R2 = rep(0,1000)
temp = 1
for(k in 1:(nrow(prediction.result.sub)/2)){
  relative_R2[temp] = as.numeric((prediction.result.sub[2*temp,"R2"]-prediction.result.sub[2*temp-1,"R2"])/prediction.result.sub[2*temp-1,"R2"])
  trait[temp] = as.character(prediction.result.sub[2*temp,"trait"])
  eth[temp] = as.character(prediction.result.sub$eth[2*temp])
  temp = temp + 1
}

relative_R2 = relative_R2[1:(temp-1)]
trait = trait[1:(temp-1)]
eth = eth[1:(temp-1)]

result = data.frame(relative_R2,trait,eth)

bintrait = c("Any CVD","Depression",
               "SBMN",
               "Migraine Diagnosis",
               "Morning Person")
contrait = c( "Heart metabolic disease burden",
             "Height")


result %>% filter(trait %in%contrait) %>% 
  group_by(eth) %>% 
  summarise(mean(relative_R2))
result %>% filter(trait %in%bintrait) %>% 
  group_by(eth) %>% 
  summarise(mean(relative_R2))

prediction.result.sub %>% filter(trait =="Any CVD")








#realtive improvement 

prediction.result.sub = prediction.result %>% 
  mutate(R2 = ifelse(trait %in% c("Heart metabolic disease burden","Height"),
                     result,
                     ifelse(result<0.5,0,qnorm(result)^2*2))) %>% 
  filter(method_vec%in%c("TDLD-SLEB","TDLD-SLEB (all ethnics)"))
trait = rep("c",1000)
eth = rep("c",1000)
relative_R2 = rep(0,1000)
temp = 1
for(k in 1:(nrow(prediction.result.sub)/2)){
  relative_R2[temp] = as.numeric((prediction.result.sub[2*temp,"R2"]-prediction.result.sub[2*temp-1,"R2"])/prediction.result.sub[2*temp-1,"R2"])
  trait[temp] = as.character(prediction.result.sub[2*temp,"trait"])
  eth[temp] = as.character(prediction.result.sub$eth[2*temp])
  temp = temp + 1
}

relative_R2 = relative_R2[1:(temp-1)]
trait = trait[1:(temp-1)]
eth = eth[1:(temp-1)]

result = data.frame(relative_R2,trait,eth)

bintrait = c("Any CVD","Depression",
             "SBMN",
             "Migraine Diagnosis",
             "Morning Person")
contrait = c( "Heart metabolic disease burden",
              "Height")


result.sub = result %>% filter(trait %in%contrait)
mean(result.sub$relative_R2)




result.sub = result %>% filter(trait %in%bintrait) 


mean(result.sub$relative_R2)





%>% 
  group_by(eth) %>% 
  summarise(mean(relative_R2))

prediction.result.sub %>% filter(trait =="Any CVD")





