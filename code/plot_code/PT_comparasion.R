#generate P + T comaprasion plot
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
for(i1 in 1:2){
  #standard LD clumping results
  load(paste0("LD.clump.result.GA_",i1,".rdata"))
  LD.clump.result <- LD.result.list[[1]]
  sample_size <- factor(rep(c("15000","45000","80000","100000"),15),
                        levels=c("15000","45000","80000","100000"))
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
  
  LD.clump.result <- cbind(LD.clump.result,sample_size,cau_vec)
  #save(LD.clump.result,file = "LD.clump.result_090420_P+T.rdata")
  #load("LD.clump.result_090420_P+T.rdata")
  
  
  LD.clump.result$eth.vec <- factor(LD.clump.result$eth.vec,
                                    #levels =c("EUR","AFR","AMR","EAS","SAS"))
                                    levels =c("AFR","AMR","EAS","EUR","SAS"))
  p <- ggplot(LD.clump.result,aes(x= sample_size,y=r2.vec,group=eth.vec))+
    geom_line(aes(color=eth.vec))+
    geom_point(aes(color=eth.vec))+
    theme_Publication()+
    ylab("R2")+
    xlab("Training sample size")+
    labs(color = "Ethnic group")+
    facet_grid(cols=vars(cau_vec))+
    scale_colour_Publication()+
    # scale_color_nejm()+
    ggtitle("C + T prediction performance for five ethnic groups")+
    theme(axis.text = element_text(face = "bold"),
          legend.text = element_text(face = "bold"))
  p
  png(file = paste0("./LD_clumping_result_summary_GA_",i1,".png"),
      width = 9.5, height = 6, res = 300,units = "in")
  print(p)
  dev.off()
}

#Underlying heritability


