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
method_vec = c("PT","BESTEUR","BESTEURLDPred2","weighted_PRS","tdld","TDLD_EB","TDLD_SLEB","TDLD_SLEBall")
method_name = c("C+T","Best EUR SNP (C+T)",
                "Best EUR SNP (LDpred2)",
                "Weighted PRS",
                "TDLD",
                "TDLD-EB",
                "CT-SLEB (two ethnics)",
                "CT-SLEB (five ethnics)")
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
    for(i in 1:5){
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
#load LDPred2 single ethic results
load("/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/R2.ldpred2.RData")
eth_group = c("european","african_american",
              "latino","east_asian","south_asian")
eth_name = c("European","African American",
             "Latino","East Asian","South Asian")
trait_name = c("Any CVD","Depression",
               "Heart metabolic disease burden",
               "Height",
               "SBMN",
               "Migraine Diagnosis",
               "Morning Person")
LDpred2 = R2.ldpred2 %>% 
  mutate(eth = case_when(
    Race=="AFR" ~ "African American",
    Race=="EUR" ~ "European",
    Race=="AMR" ~ "Latino",
    Race=="EAS" ~ "East Asian",
    Race=="SAS" ~ "South Asian"
    
  ),
  trait = case_when(
    Trait=="any_cvd" ~ "Any CVD",
    Trait=="depression" ~ "Depression",
    Trait=="heart_metabolic_disease_burden" ~ "Heart metabolic disease burden",
    Trait=="height" ~ "Height",
    Trait=="iqb.sing_back_musical_note" ~ "SBMN",
    Trait=="migraine_diagnosis" ~ "Migraine Diagnosis",
    Trait=="morning_person" ~ "Morning Person"
  ),
  result = R2,
  method_vec = "LDpred2") %>% 
  select(result,eth,trait,method_vec)
# LDpred2 = LDpred2 %>% 
#   filter(eth!="European")

prediction.result = rbind(prediction.result,LDpred2)
#load prs-csx result
prs.csx = read.csv("../PRS-CSx.csv",header = T)

prediction.result = rbind(prediction.result,prs.csx)

prediction.result = prediction.result %>% 
  filter(method_vec%in%c("TDLD-EB","TDLD","Best EUR SNP + EB coefficients (C+T)",
                         "Best EUR SNP + target coefficients (C+T)")==F)
prediction.result = prediction.result %>% 
  filter(method_vec!="PRS-CSx")
prediction.result$method_vec = factor(prediction.result$method_vec,
                                      levels = c("C+T",
                                                 "LDpred2",
                                                 "Best EUR SNP (C+T)",
                                                 "Best EUR SNP (LDpred2)",
                                                 "Weighted PRS",
                               #                  "PRS-CSx",
                                                 "CT-SLEB (two ethnics)",
                                                 "CT-SLEB (five ethnics)"))
prediction.result$Method = prediction.result$method_vec
prediction.result$index = rep("1",nrow(prediction.result))
prediction.result$eth = factor(prediction.result$eth,
                       levels = c("European","African American",
                                  "Latino","East Asian","South Asian"))

save(prediction.result,file = "/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/prediction_summary.rdata")
write.csv(prediction.result, file = "/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/prediction_summary.csv")
sigma2toauc = function(x){
  ifelse(x==0,0.50,round(pnorm(0.5*sqrt(x)),2))
}
prediction.result.table = prediction.result %>% 
  mutate(sigma2 = ifelse(result<0.5,0,qnorm(result)^2*2)) %>% 
  mutate(sigma2 = ifelse(trait%in%
                           c("Heart metabolic disease burden","Height"),NA,sigma2)) %>% 
  select(eth,trait,Method,result,sigma2)
# write.csv(prediction.result.table,file = "/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/prediction_summary.csv",row.names = F)

uvals = factor(c("C+T",
                          "LDpred2",
                          "Best EUR SNP (C+T)",
                          "Best EUR SNP (LDpred2)",
                          "Weighted PRS",
                        #  "PRS-CSx",
                 "CT-SLEB (two ethnics)",
                 "CT-SLEB (five ethnics)"),
               levels= c("C+T",
                         "LDpred2",
                         "Best EUR SNP (C+T)",
                         "Best EUR SNP (LDpred2)",
                         "Weighted PRS",
                         #"PRS-CSx",
                         "CT-SLEB (two ethnics)",
                         "CT-SLEB (five ethnics)"))

n.single = 9


single.color =  brewer.pal(n.single, "Blues")[c(4,7)]
n.EUR = 9


EUR.color = brewer.pal(n.EUR, "Greens")[c(4,7)]


n.multi = 9
#multi.color = brewer.pal(n.multi, "Oranges")[c(3,5,7,9)]
multi.color = brewer.pal(n.multi, "Oranges")[c(3,7,9)]
colour = c(single.color,EUR.color,multi.color)
col_df = tibble(
  colour = c(single.color,EUR.color,multi.color),
  method_vec = uvals,
  category = case_when(method_vec%in%c("C+T",
                                       "LDpred2") ~ "Single ethnic method",
                       method_vec%in%c("Best EUR SNP (C+T)",
                                       "Best EUR SNP (LDpred2)"
                       ) ~ "EUR PRS based method",
                       method_vec%in%c("Weighted PRS",
                                       "PRS-CSx",
                                       "CT-SLEB (two ethnics)",
                                       "CT-SLEB (five ethnics)") ~ "Multi ethnic method")
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
                      "Morning Person")) %>% 
  filter(eth!="European") %>% 
  mutate(eth = factor(eth,
                    levels = c("African American",
                                          "Latino","East Asian","South Asian")))

#create benchmark C+T or LDPred2 results for european to plot for other ethnic group
prediction.result.European.CT = prediction.result %>% 
  filter(trait %in% c("Any CVD","Depression",
                      "SBMN",
                      "Migraine Diagnosis",
                      "Morning Person")) %>% 
  filter(eth=="European"&
           method_vec=="C+T") %>% 
  select(result,trait)
prediction.result.European.LDpred2 = prediction.result %>% 
  filter(trait %in% c("Any CVD","Depression",
                      "SBMN",
                      "Migraine Diagnosis",
                      "Morning Person")) %>% 
  filter(eth=="European"&
           method_vec=="LDpred2") %>% 
  
  select(result,trait)
prediction.result.European =prediction.result %>% 
  filter(trait %in% c("Any CVD",
                      "Depression",
                      "SBMN",
                      "Migraine Diagnosis",
                      "Morning Person")) %>% 
  filter(eth=="European"&
           method_vec%in%c("C+T","LDpred2")) %>% 
  group_by(trait) %>% 
  summarise(newresult = max(result))
#replace result by better results in CT and LDpred2
 
  

# sigma2toauc = function(x){
#   ifelse(x==0,0.50,round(pnorm(0.5*sqrt(x)),2))
# }
prediction.result.sub = prediction.result.sub %>%
  mutate(sigma2 = ifelse(result<0.5,0,qnorm(result)^2*2))

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
  geom_hline(data = prediction.result.European.CT, aes(yintercept = result), linetype = "dashed",color = "red")

print(p.null)
p = plot_grid(p.null,p.leg,nrow=1,rel_widths = c(3,1))
png(filename = "../comparasion_plot/bin_comparasion.png",width=17,height = 12,units = "in",res = 300)
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


                                               
prediction.result.European = prediction.result %>% 
  filter(trait %in% c("Heart metabolic disease burden","Height")) %>% 
  filter(eth=="European") %>% 
  group_by(trait) %>% 
  summarize(newresult = max(result))
prediction.result.sub = prediction.result %>% 
  filter(trait %in% c("Heart metabolic disease burden","Height")) %>% 
  filter(eth!="European") %>% 
  mutate(eth = factor(eth,
                      levels = c("African American",
                                 "Latino","East Asian","South Asian")))

p.null <- ggplot(prediction.result.sub)+
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
  scale_fill_manual(values = colour) +
  theme(legend.position = "none")+
  geom_hline(data = prediction.result.European, aes(yintercept = newresult), linetype = "dashed",color = "red")



p = plot_grid(p.null,p.leg,nrow=1,rel_widths = c(3,1))
#theme(strip.text.y = element_blank())
#library(RColorBrewer)
png(filename = "../comparasion_plot/con_comparasion.png",width=17,height = 12,units = "in",res = 300)
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
