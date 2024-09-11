setwd("/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/LD_simulation_GA/LD_stack/")
source("../../../code/LD_simulation_large/theme_Publication.R")
library(ggplot2)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(tidyverse)
library(data.table)
load(paste0("LD.clump.result.CT.rdata"))
LD.clump.result <- LD.result.list[[1]] %>% 
  mutate(method_vec = rep("C+T")) %>% 
  filter(eth.vec == "EUR"&
           m_vec== 4)
load("LDpred2.result.021622")
LDpred2.result.sub = LDpred2.result %>% 
  filter(eth.vec=="EUR"&m_vec==4)
EUR.result = rbind(LD.clump.result,LDpred2.result.sub) %>% 
  group_by(l_vec,ga_vec) %>% 
  summarise(eur.result = max(r2.vec))

rbind(LD.clump.result,LDpred2.result.sub) %>% 
  filter(l_vec==3&ga_vec==3)
load(paste0("LD.clump.result.allethtest.EB.rdata"))
EB.result = alleth.EB.result %>% 
  mutate(method_vec = "CT-SLEB (five ancestries)")

prediction.result <- rbind(LDpred2.result.sub,
                           EB.result)


prediction.result = prediction.result %>% 
  mutate(cau_vec = case_when(
    l_vec== 1 ~ paste0("Causal SNPs Proportion = 0.01"),
    l_vec== 2 ~ paste0("Causal SNPs Proportion = 0.001"),
    l_vec== 3 ~ paste0("Causal SNPs Proportion = 5E-04")
  ),
  sample_size = case_when(
    m_vec == 1~ "15000",
    m_vec == 2~ "45000",
    m_vec == 3~ "80000",
    m_vec == 4~ "100000"
  )) %>% 
  mutate(cau_vec = factor(cau_vec,
                          levels = c("Causal SNPs Proportion = 0.01",
                                     "Causal SNPs Proportion = 0.001",
                                     "Causal SNPs Proportion = 5E-04")),
         sample_size = factor(sample_size,
                              levels = c("15000","45000","80000","100000")))
i1 = 1
  result.sub.list = list()
  l = 2
    result.sub = prediction.result %>% 
      filter(l_vec==l&
               ga_vec ==i1&
               eth.vec!="EUR")
    result.sub.EUR = prediction.result %>% 
      filter(l_vec==l&
               ga_vec ==i1&
               eth.vec=="EUR")
    result.sub$relative.r2 = result.sub$r2.vec/result.sub.EUR$r2.vec
    result.sub.list[[l]] = result.sub
  
  result.sub = rbindlist(result.sub.list)
  p <- ggplot(result.sub,aes(x= sample_size,y=relative.r2,group=eth.vec))+
    geom_line(aes(color=eth.vec),size = 1)+
    geom_point(aes(color=eth.vec))+
    ylab(expression(bold(paste(R^2, "of CT-SLEB in the target population / ",R^2," of EUR PRS in EUR"))))+
    theme_Publication()+
    xlab("Training sample size")+
    labs(color = "Ancestry group")+
    facet_grid(cols=vars(cau_vec))+
    #scale_color_brewer(palette = "Paired")+
    scale_colour_Publication()+
    ggtitle(paste0("Relative performances of CT-SLEBover the LDpred2 EUR PRS in EUR"))+
    geom_hline(yintercept=1, linetype="dashed", 
               color = "red")+
   theme(legend.text=element_text(size=12),
          axis.text = element_text(size=12),
          axis.title.y = element_text(angle=90,vjust =2,size = 12),
         panel.grid.major = element_line(colour="#f0f0f0"),
         panel.grid.minor = element_line(colour="#f0f0f0"))+
    theme(plot.title = element_text(face = "bold",
                                    size = rel(1.0), hjust = 0.5))
  print(p)
  if(i1 ==1){
    p1 = p
  }else if(i1==2){
    p2 = p
  }
  
  png(file = paste0("/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/presentation_plot/presentation_fix_gap_summary_GA_",i1,".png"),
      width = 10, height = 7, res = 300,units = "in")
  print(p)
  dev.off()




  
  i1 = 2
  result.sub.list = list()
  l = 2
  result.sub = prediction.result %>% 
    filter(l_vec==l&
             ga_vec ==i1&
             eth.vec!="EUR")
  result.sub.EUR = prediction.result %>% 
    filter(l_vec==l&
             ga_vec ==i1&
             eth.vec=="EUR")
  result.sub$relative.r2 = result.sub$r2.vec/result.sub.EUR$r2.vec
  result.sub.list[[l]] = result.sub
  
  result.sub = rbindlist(result.sub.list)
  p <- ggplot(result.sub,aes(x= sample_size,y=relative.r2,group=eth.vec))+
    geom_line(aes(color=eth.vec),size = 1)+
    geom_point(aes(color=eth.vec))+
    ylab(expression(bold(paste(R^2, "of CT-SLEB in the target population / ",R^2," of EUR PRS in EUR"))))+
    theme_Publication()+
    xlab("Training sample size")+
    labs(color = "Ancestry group")+
    facet_grid(cols=vars(cau_vec))+
    #scale_color_brewer(palette = "Paired")+
    scale_colour_Publication()+
    ggtitle(paste0("Relative performances of CT-SLEBover the LDpred2 EUR PRS in EUR"))+
    geom_hline(yintercept=1, linetype="dashed", 
               color = "red")+
    theme(legend.text=element_text(size=12),
          axis.text = element_text(size=12),
          axis.title.y = element_text(angle=90,vjust =2,size = 12),
          panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.minor = element_line(colour="#f0f0f0"))+
    theme(plot.title = element_text(face = "bold",
                                    size = rel(1.0), hjust = 0.5))
  print(p)
  if(i1 ==1){
    p1 = p
  }else if(i1==2){
    p2 = p
  }
  
  png(file = paste0("/Users/zhangh24/Library/CloudStorage/Box-Box/multi_ethnic/result/presentation_plot/presentation_fix_gap_summary_GA_",i1,".png"),
      width = 10, height = 7, res = 300,units = "in")
  print(p)
  dev.off()
  
  