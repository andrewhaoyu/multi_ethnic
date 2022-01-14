setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA/LD_stack/")
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
load(paste0("LDpred2.result.rdata"))
LDpred2.result.sub = LDpred2.result %>% 
  filter(eth.vec=="EUR"&m_vec==4)
EUR.result = rbind(LD.clump.result,LDpred2.result.sub) %>% 
  group_by(l_vec,ga_vec) %>% 
  summarise(eur.result = max(r2.vec))

rbind(LD.clump.result,LDpred2.result.sub) %>% 
  filter(l_vec==3&ga_vec==3)
load(paste0("LD.clump.result.EB.rdata"))
EB.result = EB.result %>% 
  filter(method_vec=="TDLD-SLEB")
prediction.result <- rbind(LD.clump.result,
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
for(i1 in 1:5){
  result.sub.list = list()
  for(l in 1:3){
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
  }
   result.sub = rbindlist(result.sub.list)
  p <- ggplot(result.sub,aes(x= sample_size,y=relative.r2,group=eth.vec))+
    geom_line(aes(color=eth.vec),size = 1)+
    geom_point(aes(color=eth.vec))+
    ylab(expression(bold(paste("Relative ", R^2))))+
    theme_Publication()+
    xlab("Training sample size")+
    labs(color = "Ethnic group")+
    facet_grid(cols=vars(cau_vec))+
    #scale_color_brewer(palette = "Paired")+
    scale_colour_Publication()+
    ggtitle(paste0("Relative performance of CT-SLEB over single ethnic European PRS"))+
    geom_hline(yintercept=1, linetype="dashed", 
               color = "red")+
    theme(legend.text=element_text(size=12))
  print(p)
  if(i1 ==1){
    p1 = p
  }else if(i1==2){
    p2 = p
  }
 
  png(file = paste0("./fix_gap_summary_GA_",i1,".png"),
      width = 12.6, height = 7, res = 300,units = "in")
  print(p)
  dev.off()
  }







LD.clump.result = LD.clump.result %>% 
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
                              levels = c("15000","45000","80000","100000")),
  ) %>% 
  mutate(ga_arc = case_when(ga_vec==1 ~"Fixed common SNP heritability",
                            ga_vec==2 ~"Strong negative selection",
                            ga_vec==3 ~"Less correlation",
                            ga_vec==4 ~"No negative selection",
                            ga_vec==5 ~"Mild negative selection",
  ))
LD.clump.result.plot = LD.clump.result %>% 
  filter(ga_vec==1)
p <- ggplot(LD.clump.result.plot,aes(x= sample_size,y=r2.vec,group=eth.vec))+
  geom_line(aes(color=eth.vec))+
  geom_point(aes(color=eth.vec))+
  theme_Publication()+
  ylab("R2")+
  xlab("Training sample size")+
  labs(color = "Ethnic group")+
  facet_grid(cols=vars(cau_vec))+
  #scale_color_brewer(palette = "Paired")+
  scale_color_npg()+
  ggtitle("Prediction performance across ethnic groups using P+T")


p <- ggplot(LD.clump.result.plot,aes(x= sample_size,y=r2.vec,group=eth.vec))+
  geom_line(aes(color=eth.vec))+
  geom_point(aes(color=eth.vec))+
  theme_Publication()+
  ylab("R2")+
  xlab("Training sample size")+
  labs(color = "Ethnic group")+
  facet_grid(cols=vars(cau_vec))+
  #scale_color_brewer(palette = "Paired")+
  scale_color_npg()+
  ggtitle("Prediction performance across ethnic groups using P+T")

setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA")
# library(ggplot2)
# library(ggsci)
#standard heritability
herita.table
# sample_size <- factor(rep(c("15000","45000","80000","100000"),6),
#                                            levels=c("15000","45000","80000","100000"))
cau_vec <- as.character(LD.clump.result$l_vec)
csp <- c(0.01,0.001,0.0005)
for(l in 1:3){
  idx <- which(LD.clump.result$l_vec==l)
  cau_vec[idx] <- paste0("Causal SNPs Proportion = ",csp[l])
}
cau_vec <- factor(cau_vec,
                  levels = paste0("Causal SNPs Proportion = ",csp))

herita.table <- cbind(herita.table,sample_size,cau_vec)
# save(LD.clump.result,file = "LD.clump.result_090420_P+T.rdata")
# load("LD.clump.result_090420_P+T.rdata")


herita.table$eth_vec <- factor(herita.table$eth_vec,
                               #levels =c("EUR","AFR","AMR","EAS","SAS"))
                               levels =c("EUR","AFR","AMR","EAS","SAS"))
p <- ggplot(herita.table,aes(x= cau_vec,y=herit_vec))+
  geom_bar(aes(fill=eth_vec),
           stat="identity",
           position = position_dodge())+
  theme_Publication()+
  ylab("R2")+
  xlab("Training sample size")+
  labs(color = "Ethnic group")+
  #facet_grid(cols=vars(cau_vec))+
  scale_fill_nejm()+
  ggtitle("Underlying heritability for different ethnic groups")





custom.col <- c("#FFDB6D", "#C4961A", "#F4EDCA", 
                "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
