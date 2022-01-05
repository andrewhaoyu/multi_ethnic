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
load(paste0("weightedprs.result.rdata"))
load(paste0("LD.clump.result.2DLD.rdata"))
load(paste0("LD.clump.result.EB.rdata"))
load(paste0("LD.clump.result.alleth.EB.rdata"))
load(paste0("LDpred2.result.rdata"))
load(paste0("LDpredEUR.result.rdata"))
load(paste0("prscsx.result.rdata"))
LD.clump.result <- LD.result.list[[1]] %>% 
  mutate(method_vec = rep("C+T"))

TDLD.result = TDLD.result %>% 
  filter(method_vec=="TDLD")

alleth.EB.result = alleth.EB.result %>% 
  filter(method_vec=="TDLD-SLEB (all ethnics)")
weightedprs.result = weightedprs.result %>% 
  mutate(method_vec = "Weighted PRS")
EB.result = EB.result %>% 
  mutate(method_vec=ifelse(method_vec=="TDLD-SLEB","CT-SLEB (two ethnics)",method_vec))
alleth.EB.result = alleth.EB.result %>% 
  mutate(method_vec=ifelse(method_vec=="TDLD-SLEB (all ethnics)","CT-SLEB (five ethnics)",method_vec))
prediction.result <- rbind(LD.clump.result,
                           SCT.clump.result,
                           LDpred2.result,
                           eursnp.result,
                           LDpredEUR.result,
                           weightedprs.result,
                           prscsx.result,
                           TDLD.result,
                           EB.result,
                           alleth.EB.result)


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
                              levels = c("15000","45000","80000","100000")),
         method_vec = factor(method_vec,
                             levels = c("C+T",
                                        "SCT",
                                        "LDpred2",
                                        "Best EUR SNP (C+T)",
                                        "Best EUR SNP + target coefficients (C+T)",
                                        "Best EUR SNP + EB coefficients (C+T)",
                                        "Best EUR PRS (LDpred2)",
                                        "Weighted PRS",
                                        "PRS-CSx",
                                        "TDLD",
                                        "TDLD-EB",
                                        "CT-SLEB (two ethnics)",
                                        "CT-SLEB (five ethnics)"
                             ))) %>% 
  mutate(ga_arc = case_when(ga_vec==1 ~"Fixed common SNP heritability with strong negative selection",
                            ga_vec==2 ~"Fixed whole genome heritability with strong negative selection",
                            ga_vec==3 ~"Fixed common SNP heritability with strong negative selection (r=0.6)",
                            ga_vec==4 ~"Fixed common SNP heritability with no negative selection",
                            ga_vec==5 ~"Fixed common SNP heritability with mild negative selection",
  ))

prediction.result = prediction.result %>% 
  filter(method_vec%in%
           c("C+T",
             "LDpred2",
             "Best EUR SNP (C+T)",
             "Best EUR PRS (LDpred2)",
             "Weighted PRS",
             "PRS-CSx",
             "CT-SLEB (two ethnics)",
             "CT-SLEB (five ethnics)"
           ))

uvals = unique(prediction.result$method_vec)

n.single = 9


single.color =  brewer.pal(n.single, "Blues")[c(4,7)]
n.EUR = 9


EUR.color = brewer.pal(n.EUR, "Greens")[c(4,7)]


n.multi = 9
multi.color = brewer.pal(n.multi, "Oranges")[c(3,5,7,9)]
colour = c(single.color,EUR.color,multi.color)
col_df = tibble(
  colour = c(single.color,EUR.color,multi.color),
  method_vec = uvals,
  category = case_when(method_vec%in%c("C+T",
                                       "LDpred2") ~ "Single ethnic method",
                       method_vec%in%c("Best EUR SNP (C+T)",
                                       "Best EUR PRS (LDpred2)"
                       ) ~ "EUR PRS based method",
                       method_vec%in%c("Weighted PRS",
                                       "PRS-CSx",
                                       "CT-SLEB (two ethnics)",
                                       "CT-SLEB (five ethnics)") ~ "Multi ethnic method")
) %>% 
  mutate(category = factor(category,levels = c("Single ethnic method",
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
    prediction.result.sub %>% 
      filter(category %in% filler),
    aes(x= sample_size,y=r2.vec,
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

    m =1
    i1 = 1
    prediction.result.sub <- prediction.result %>% 
      filter(ga_vec==i1&
               eth.vec!="EUR"&
               m_vec ==m&
               l_vec==3&
               eth.vec=="EAS")
    title = prediction.result.sub$ga_arc[1]
    legs = lapply(sort(unique(col_df$category)), run_plot)
    
    legs = lapply(legs, getLegend)
    p.leg = plot_grid(NULL,NULL,legs[[1]],legs[[2]], legs[[3]],NULL,align="v",ncol=1,rel_heights=c(1,1,0.8,0.9,1.2,1.2))
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
        width = 10, height = 8, res = 300,units = "in")
    print(p.null)
    dev.off()
    
    prediction.result.sub <- prediction.result %>% 
      filter(ga_vec==i1&
               eth.vec!="EUR"&
               m_vec ==3&
               l_vec==3&
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
    png(file = paste0("../../presentation_plot/single_ethnic_single_cau_compare_methods_large.png"),
        width = 14, height = 8, res = 300,units = "in")
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
        width = 14, height = 8, res = 300,units = "in")
    print(p)
    dev.off()
    
    
    