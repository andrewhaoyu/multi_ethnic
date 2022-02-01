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
method_vec = c("TDLD_SLEB","TDLD_SLEBall")
method_name = c(
                "old CT-SLEB (two ancestries)",
                "old CT-SLEB (five ancestries)")
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
      if(i==1&i1==1){
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
      }else if(i==1&i1!=1){
        #do nothing since European only have single ethnic method 
      }else {
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
        }else if(method_vec[i1]=="BESTEURLDPred2"){
          result = read.csv(paste0("./ldpred2_EUR/",eth_group[i],"_",trait[l],"_ldpred2_baseline_testing"))
          plot.data = data.frame(result = result,
                                 eth = eth_name[i],
                                 trait = trait_name[l],
                                 method_vec = method_name[i1])
          colnames(plot.data)[1] = "result"
          plot.data.list[[temp]] = plot.data
          temp = temp+1
        }
        else{
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
}
prediction.result = rbindlist(plot.data.list)
tdld = read.csv("../TDLD_EB.csv",header = T)

prediction.result = rbind(prediction.result,tdld)

prediction.result$method_vec = factor(prediction.result$method_vec,
                                      levels = c("old CT-SLEB (two ancestries)",
                                                 "CT-SLEB (two ancestries)",
                                                 "old CT-SLEB (five ancestries)",
                                                 "CT-SLEB (five ancestries)"))
prediction.result$Method = prediction.result$method_vec
prediction.result$index = rep("1",nrow(prediction.result))
prediction.result$eth = factor(prediction.result$eth,
                               levels = c("European","African American",
                                          "Latino","East Asian","South Asian"))



sigma2toauc = function(x){
  ifelse(x==0,0.50,round(pnorm(0.5*sqrt(x)),2))
}
prediction.result.table = prediction.result %>% 
  mutate(sigma2 = ifelse(result<0.5,0,qnorm(result)^2*2)) %>% 
  mutate(sigma2 = ifelse(trait%in%
                           c("Heart metabolic disease burden","Height"),NA,sigma2)) %>% 
  select(eth,trait,Method,result,sigma2)


prediction.result.sub = prediction.result %>% 
  filter(trait %in% c("Any CVD","Depression",
                      "SBMN",
                      "Migraine Diagnosis",
                      "Morning Person")) %>% 
  filter(eth!="European") %>% 
  mutate(eth = factor(eth,
                      levels = c("African American",
                                 "Latino","East Asian","South Asian")))


p <- ggplot(prediction.result.sub)+
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
  scale_fill_Publication()
  #scale_fill_manual(values = colour) +
  #theme(legend.position = "none")+
 # geom_hline(data = prediction.result.European.CT, aes(yintercept = result), linetype = "dashed",color = "red")



png(filename = "../comparasion_plot/CT_SLEB_old_new_bin.png",width=17,height = 12,units = "in",res = 300)
print(p)
dev.off()

# #convert AUC to liability scale variance
# #AUC = pnorm(sqrt(sigma^2/2))
# #sigma^2 = qnorm(AUC)^2*2
# prediction.result.sub = prediction.result.sub %>% 
#   mutate(sigma2 = ifelse(result<0.5,0,qnorm(result)^2*2))
# 
# p.null <- ggplot(prediction.result.sub)+
#   geom_bar(aes(x = index,y = sigma2,fill=method_vec),
#            position = position_dodge(),
#            stat = "identity")+
#   theme_Publication()+
#   ylab("Liability scale variance")+
#   facet_grid(vars(trait),vars(eth))+
#   theme_Publication()+
#   #coord_cartesian(ylim = c(0.47, 0.65)) +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())+
#   scale_fill_manual(values = colour) +
#   theme(legend.position = "none")
# 
# p = plot_grid(p.null,p.leg,nrow=1,rel_widths = c(3,1))
# png(filename = "../comparasion_plot/bin_comparasion_r2.png",width=17,height = 11,units = "in",res = 300)
# print(p)
# dev.off()



prediction.result.sub = prediction.result %>% 
  filter(trait %in% c("Heart metabolic disease burden","Height")) %>% 
  filter(eth!="European") %>% 
  mutate(eth = factor(eth,
                      levels = c("African American",
                                 "Latino","East Asian","South Asian")))

p <- ggplot(prediction.result.sub)+
  geom_bar(aes(x = index,y = result,fill=method_vec),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  ylab(expression(bold(R^2)))+
  facet_grid(vars(trait),vars(eth),scales = "free")+
  theme_Publication()+
  #coord_cartesian(ylim = c(0.47, 0.65)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_Publication()



#theme(strip.text.y = element_blank())
#library(RColorBrewer)
png(filename = "../comparasion_plot/CT_SLEB_old_new_con.png",width=17,height = 12,units = "in",res = 300)
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

result %>% group_by(eth) %>% 
  summarise(mean(relative_R2))
prediction.result.sub %>% filter(trait =="Any CVD")
