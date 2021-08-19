#plot the bar plot for different chips P + T comparasion
setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA")
source("../../code/LD_simulation_large/theme_Publication.R")
library(ggplot2)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
#library(RColorBrewer)
colourCount = 9
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#
for(i1 in 2:2){
  
  #standard LD clumping results
  load(paste0("LD.clump.result.GA_",i1,".rdata"))
  LD.clump.result.1kg <- LD.result.list[[1]]
  chips = rep("1kg",nrow(LD.clump.result.1kg))
  LD.clump.result.1kg$chips = chips
  load(paste0("LD.clump.result.GA_",i1,"_mega.rdata"))
  LD.clump.result.mega <- LD.result.list[[1]]
  chips = rep("mega",nrow(LD.clump.result.mega))
  LD.clump.result.mega$chips = chips
  load(paste0("LD.clump.result.GA_",i1,"_hm.rdata"))
  LD.clump.result.hm <- LD.result.list[[1]]
  chips = rep("hm",nrow(LD.clump.result.mega))
  LD.clump.result.hm$chips = chips
  
  LD.clump.result = rbind(LD.clump.result.1kg,LD.clump.result.mega,LD.clump.result.hm)
  
  #save(LD.clump.result,file = "LD.clump.result_090420_P+T.rdata")
  #load("LD.clump.result_090420_P+T.rdata")
  sample_size =  as.character(LD.clump.result$method_vec)
  ssp = as.character(c("15000","45000","80000","100000"))
  cau_vec <- as.character(LD.clump.result$l_vec)
  csp <- c(0.01,0.001,0.0005)
  for(l in 1:3){
    idx <- which(LD.clump.result$l_vec==l)
    cau_vec[idx] <- paste0("Causal SNPs Proportion = ",csp[l])
  }
  cau_vec = factor(cau_vec,levels = paste0("Causal SNPs Proportion = ",csp))
  for(m in 1:4){
    idx <- which(LD.clump.result$m_vec==m)
    sample_size[idx] <- ssp[m]
  }
  sample_size = factor(sample_size,
                       levels = c("15000","45000","80000","100000"))
  LD.clump.result <- cbind(LD.clump.result,cau_vec,sample_size)
  
  
  LD.clump.result$eth.vec <- factor(LD.clump.result$eth.vec,
                                    #levels =c("EUR","AFR","AMR","EAS","SAS"))
                                    levels =c("AFR","AMR","EAS","EUR","SAS")) 
    LD.clump.result = LD.clump.result %>%  
    mutate(
      Arrays = case_when(chips=="1kg" ~ "1KG",
                              chips=="hm" ~ "HM3",
                              chips=="mega"~"HM3 + MEGA"))
  
  for(m in 1:4){
    idx <- which(LD.clump.result$m_vec ==m)
    LD.clump.result.sub = LD.clump.result[idx,]
    p <- ggplot(LD.clump.result.sub,aes(x= sample_size,y=r2.vec,group=chips))+
      geom_bar(aes(fill=Arrays),
               stat="identity",
               position = position_dodge())+
      #geom_point(aes(color=method_vec))+
      theme_Publication()+
      ylab("R2")+
      xlab("Sample Size")+
      labs(fill = "Method")+
      facet_grid(vars(cau_vec),vars(eth.vec))+
      scale_fill_nejm()+
      theme(axis.text = element_text(size = rel(0.9)),
            legend.text = element_text(size = rel(0.9)))+
      ggtitle("Prediction performance comparasion (C + T)")
    p
    
    png(file = paste0("./LD_clumping_compare_GA_",i1,"_size_",m,".png"),
        width = 10, height = 8, res = 300,units = "in")
    print(p)
    dev.off()
    
  }
  
}









