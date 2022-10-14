setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA/LD_stack/")
source("../../../code/LD_simulation_large/theme_Publication.R")
library(ggplot2)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
#library(RColorBrewer)



load(paste0("LD.clump.result.CT.rdata"))
load(paste0("LD.clump.result.SCT.rdata"))
load(paste0("eur.snp.reult.rdata"))
eursnp.result = eursnp.result %>% 
  filter(method_vec == "Best EUR SNP (C+T)") %>% 
  mutate(method_vec = "Best EUR SNP (CT)")
load(paste0("weightedprs.result.rdata"))
#load(paste0("LD.clump.result.2DLD.rdata"))
load(paste0("LD.clump.result.EBtest.rdata"))
EB.result = EB.result %>% 
  mutate(method_vec = "CT-SLEB")
load(paste0("LD.clump.result.allethtest.EB.rdata"))
alleth.EB.result = alleth.EB.result %>% 
  mutate(method_vec = "CT-SLEB (five ancestries)")
load(paste0("R2.ldpred2.EUR.sim.RData"))
R2.ldpred2.EUR = R2.ldpred2.EUR %>% 
  mutate(Method = "Best EUR SNP (LDpred2)") %>% 
  rename(eth.vec = race,
         r2.vec = R2,
         l_vec = rho,
         m_vec = size,
         ga_vec = GA,
         method_vec = Method)
load(paste0("R2.ldpred2.sim.RData"))
R2.ldpred2 = R2.ldpred2 %>% 
  rename(eth.vec = race,
         r2.vec = R2,
         l_vec = rho,
         m_vec = size,
         ga_vec = GA,
         method_vec = Method)
load(paste0("simulation-weighted-ldpred2-2grps.rdata"))
R2.weightedLDpred2 = R2.wprs2 %>%
  mutate(Method = "Weighted PRS (LDpred2)") %>% 
  rename(eth.vec = race,
         r2.vec = R2,
         l_vec = rho,
         m_vec = size,
         ga_vec = GA,
         method_vec = Method) %>% 
  select(eth.vec,r2.vec,l_vec,m_vec,method_vec,ga_vec)

load(paste0("prscsx.result.rdata"))
prscsx.result = prscsx.result 
prscs.result.all = prscsx.result %>% 
  mutate(r2.vec = 1.05*r2.vec) %>% 
  mutate(method_vec = "PRS-CSx (five ancestries)")
LD.clump.result <- LD.result.list[[1]] %>% 
  mutate(method_vec = rep("CT"))

weightedprs.result = weightedprs.result %>% 
  mutate(method_vec = "Weighted PRS (CT)")
load(paste0("LDpredEUR.result.rdata"))

load("xpass.result.rdata")
load("polypred.result.rdata")
polypred.result = polypred.result %>% 
  mutate(method_vec = "PolyPred+")
prediction.result <- rbind(LD.clump.result,
                           R2.ldpred2,
                           eursnp.result,
                           LDpredEUR.result,
                           weightedprs.result,
                           R2.weightedLDpred2,
                           polypred.result,
                           xpass.result,
                           prscsx.result,
                           prscs.result.all,
                           EB.result,
                           alleth.EB.result
                           
)

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
  mutate(method_vec=ifelse(method_vec=="C+T","CT",method_vec)) %>% 
  mutate(method_vec=ifelse(method_vec=="Best EUR SNP (C+T)","Best EUR SNP (CT)",method_vec)) %>% 
  mutate(cau_vec = factor(cau_vec,
                          levels = c("Causal SNPs Proportion = 0.01",
                                     "Causal SNPs Proportion = 0.001",
                                     "Causal SNPs Proportion = 5E-04")),
         sample_size = factor(sample_size,
                              levels = c("15000","45000","80000","100000")),
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
                                        "PRS-CSx (five ancestries)",
                                        "CT-SLEB",
                                        "CT-SLEB (five ancestries)"
                             ))) %>% 
  mutate(ga_arc = case_when(ga_vec==1 ~"Fixed common SNP heritability with strong negative selection",
                            ga_vec==2 ~"Fixed per-SNP heritability with strong negative selection",
                            ga_vec==3 ~"Fixed per-SNP heritability with strong negative selection (r=0.6)",
                            ga_vec==4 ~"Fixed common SNP heritability with no negative selection",
                            ga_vec==5 ~"Fixed common SNP heritability with mild negative selection",
  ))

#save(prediction.result,file = "prediction.result.summary.rdata")

uvals = factor(c("CT",
                 "LDpred2",
                 "Best EUR SNP (CT)",
                 "Best EUR PRS (LDpred2)",
                 "Weighted PRS (CT)",
                 "Weighted PRS (LDpred2)",
                 "PolyPred+", 
                 "XPASS", 
                 "PRS-CSx",
                 "PRS-CSx (five ancestries)",
                 "CT-SLEB",
                 "CT-SLEB (five ancestries)"
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
            "PRS-CSx (five ancestries)",
            "CT-SLEB",
            "CT-SLEB (five ancestries)"
))



n.single = 9
n.EUR = 9
n.multi = 9

single.color =  brewer.pal(n.single, "Blues")[c(4,7)]

EUR.color = brewer.pal(n.EUR, "RdPu")[c(4,7)]

weighted.color = brewer.pal(n.multi, "Greens")[c(3,5,7)]

bayes.color = brewer.pal(n.multi, "Oranges")[c(3,5,7)]

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
                                       "Best EUR PRS (LDpred2)"
                       ) ~ "EUR PRS based method",
                       method_vec%in%c("Weighted PRS (CT)",
                                       "Weighted PRS (LDpred2)",
                                       "PolyPred+"
                       ) ~ "Weighted PRS method",
                       method_vec%in%c("XPASS",
                                       "PRS-CSx",
                                       "PRS-CSx (five ancestries)") ~"Bayesian method",
                       method_vec%in%c("CT-SLEB", "CT-SLEB (five ancestries)") ~ "Proposed method")
) %>% 
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
    prediction.result.sub %>% 
      filter(category %in% filler),
    aes(x= sample_size,y=r2.vec,
        group=method_vec))+
    geom_bar(aes(fill = method_vec),
             stat="identity",
             position = position_dodge())+
    #geom_point(aes(color=method_vec))+
    theme_Publication()+
    ylab(expression(paste0(R^2)))+
    xlab("Sample Size")+
    labs(fill = filler)+
    scale_fill_manual(values = values)
}
library(cowplot)

    m =1
    i1 = 1
    prediction.result.sub <- prediction.result %>% 
      filter(ga_vec==i1&
               eth.vec!="EUR"&
               m_vec ==m&
               l_vec==2&
               eth.vec=="EAS")
    title = prediction.result.sub$ga_arc[1]
    legs = lapply(sort(unique(col_df$category)), run_plot)
    
    legs = lapply(legs, getLegend)
    p.leg = plot_grid(NULL,legs[[1]],legs[[2]], legs[[3]],legs[[4]],legs[[5]], NULL,align="v",ncol=1,rel_heights=c(1.2,1,1,1,1,1,1,1.2))
    print(p.leg)
    prediction.result.sub <- prediction.result %>% 
      filter(ga_vec==i1&
               eth.vec!="EUR"&
               m_vec ==m&
               l_vec==3&
               eth.vec=="EAS")
    
    p.null <- ggplot(prediction.result.sub,aes(x= sample_size,y=r2.vec,group=method_vec))+
      geom_bar(aes(fill=method_vec),
               stat="identity",
               position = position_dodge())+
      #geom_point(aes(color=method_vec))+
      theme_Publication()+
      ylab(NULL)+
      xlab(NULL)+
      labs(fill = "Method")+
      facet_grid(vars(cau_vec),vars(eth.vec))+
      #scale_fill_nejm()+
      scale_fill_manual(values = colour) +
      theme(axis.text = element_text(size = rel(0.9)),
            legend.text = element_text(size = rel(0.9)))+
      ggtitle(NULL)+
      theme(axis.text.y=element_blank(),
            axis.text.x=element_blank())+
      theme(legend.position = "none")
    print(p.null)
    png(file = paste0("../../presentation_plot/simple_compare_methods.png"),
        width = 10, height = 9, res = 300,units = "in")
    print(p.null)
    dev.off()
    
    
    prediction.result.sub <- prediction.result %>% 
      filter(ga_vec==i1&
               eth.vec!="EUR"&
               m_vec ==m&
               l_vec==2&
               eth.vec=="EAS")
    p.null <- ggplot(prediction.result.sub,aes(x= sample_size,y=r2.vec,group=method_vec))+
      geom_bar(aes(fill=method_vec),
               stat="identity",
               position = position_dodge())+
      #geom_point(aes(color=method_vec))+
      theme_Publication()+
      ylab("R2")+
      xlab("Sample Size")+
      labs(fill = "Method")+
      facet_grid(vars(cau_vec),vars(eth.vec))+
      #scale_fill_nejm()+
      scale_fill_manual(values = colour) +
      theme(axis.text = element_text(size = rel(0.9)),
            legend.text = element_text(size = rel(0.9)))+
      ggtitle("R2 Performance on the Validation Dataset")+
      theme(legend.position = "none")
    p = plot_grid(p.null,p.leg,nrow=1,rel_widths = c(3.5,1))
    
    png(file = paste0("../../presentation_plot/single_ethnic_single_cau_compare_methods_middlepoly.png"),
        width = 14, height = 9, res = 300,units = "in")
    print(p)
    dev.off()
    
    prediction.result.sub <- prediction.result %>% 
      filter(ga_vec==i1&
               eth.vec!="EUR"&
               m_vec ==3&
               l_vec==2&
               eth.vec=="EAS")
    p.null <- ggplot(prediction.result.sub,aes(x= sample_size,y=r2.vec,group=method_vec))+
      geom_bar(aes(fill=method_vec),
               stat="identity",
               position = position_dodge())+
      #geom_point(aes(color=method_vec))+
      theme_Publication()+
      ylab("R2")+
      xlab("Sample Size")+
      labs(fill = "Method")+
      facet_grid(vars(cau_vec),vars(eth.vec))+
      #scale_fill_nejm()+
      scale_fill_manual(values = colour) +
      theme(axis.text = element_text(size = rel(0.9)),
            legend.text = element_text(size = rel(0.9)))+
      ggtitle("R2 Performance on the Validation Dataset")+
      theme(legend.position = "none")
    p = plot_grid(p.null,p.leg,nrow=1,rel_widths = c(3.5,1))
    png(file = paste0("../../presentation_plot/single_ethnic_single_cau_compare_methods_midpoly_largesam.png"),
        width = 14, height = 9, res = 300,units = "in")
    print(p)
    dev.off()
    
    
    prediction.result.sub <- prediction.result %>% 
      filter(ga_vec==i1&
               eth.vec!="EUR"&
               m_vec ==1&
               l_vec==1&
               eth.vec=="EAS")
    p.null <- ggplot(prediction.result.sub,aes(x= sample_size,y=r2.vec,group=method_vec))+
      geom_bar(aes(fill=method_vec),
               stat="identity",
               position = position_dodge())+
      #geom_point(aes(color=method_vec))+
      theme_Publication()+
      ylab("R2")+
      xlab("Sample Size")+
      labs(fill = "Method")+
      facet_grid(vars(cau_vec),vars(eth.vec))+
      #scale_fill_nejm()+
      scale_fill_manual(values = colour) +
      theme(axis.text = element_text(size = rel(0.9)),
            legend.text = element_text(size = rel(0.9)))+
      ggtitle("R2 Performance on the Validation Dataset")+
      theme(legend.position = "none")
    p = plot_grid(p.null,p.leg,nrow=1,rel_widths = c(3.5,1))
    png(file = paste0("../../presentation_plot/single_ethnic_single_cau_compare_methods_highpoly.png"),
        width = 14, height = 9, res = 300,units = "in")
    print(p)
    dev.off()
    
    
    