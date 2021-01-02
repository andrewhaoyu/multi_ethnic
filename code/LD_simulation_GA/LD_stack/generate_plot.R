setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA")
source("../../code/LD_simulation_large/theme_Publication.R")
library(ggplot2)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
#library(RColorBrewer)
colourCount = 12
getPalette = colorRampPalette(brewer.pal(9, "Paired"))







for(i1 in 1:2){
  load(paste0("LD.clump.result.GA_",i1,".rdata"))
  LD.clump.result <- LD.result.list[[1]]
  LD.clump.result <- LD.clump.result[LD.clump.result$eth.vec!="EUR",]
  method_vec <- rep("P+T",nrow(LD.clump.result))
  LD.clump.result$method_vec = method_vec
  LD.clump.result.PT <- LD.clump.result
  #Best EUR result
  load(paste0("best_eur_snp_result_GA_",i1,".rdata"))
  method <- c("eurcoef","tarcoef","eb")
  method_nameupdate <- c("Best EUR PRS","Best EUR SNP + target coefficients","Best EUR SNP + EB")
  for(q in 1:length(method)){
    idx <- which(LD.clump.result$method_vec==method[q])
    LD.clump.result$method_vec[idx] <-   method_nameupdate[q]
  }
  LD.clump.result.EUR <- LD.clump.result
  #
  
  LD.clump.result <- rbind(LD.clump.result.EUR,LD.clump.result.PT)
  
  #"2DLD","2DLD-EB"))
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
  
  load("R2.simulation.hapmap3.RData")
  res$cau_vec = as.character(res$cau_vec)
  idx <- which(res$cau_vec=="Causal SNPs Proportion = 5e-4")
  res$cau_vec[idx] = "Causal SNPs Proportion = 5e-04"
  
  res$sample_size = as.character(res$sample_size)
  
  idx <- which(res$sample_size=="1e+05")
  res$sample_size[idx] = "100000"
  res = res %>% 
    mutate(m_vec = case_when(sample_size%in%"15000"~1,
                             sample_size%in%"45000"~2,
                             sample_size%in%"80000"~3,
                             sample_size%in%"100000"~4))
  
  res$sample_size = factor(res$sample_size,
                           levels = c("15000","45000","80000","100000"))
  
  
  LD.result.LDpred <- res %>%
    rename(eth.vec = race,
           r2.vec = R2,
           method_vec = Method) %>% 
    filter(GA==i1) %>% 
    select(eth.vec,r2.vec,cau_vec,sample_size,method_vec,m_vec)
  LD.clump.result = LD.clump.result %>% 
    select(colnames(LD.result.LDpred))
  LD.clump.result.plot = rbind(LD.clump.result,LD.result.LDpred)
  LD.clump.result.sub = LD.clump.result.plot %>% 
    filter(m_vec==1)
  save(LD.clump.result.plot,file = paste0("LD.clump.pred.result_GA_",i1,".rdata"))
  load(paste0("./LD_stack/LD.clump.result.GA_",i1,".rdata"))
  LD.result.stack = LD.result.list[[1]]
  load(paste0("./LD_stack/LD.clump.result.two_way_GA_",i1,".rdata"))
  LD.result.2D.stack = LD.result.list[[1]]
  load(paste0("./LD_stack/LD.clump.result.eb_GA_",i1,".rdata"))
  LD.result.eb.stack = LD.result.list[[1]]
  load(paste0("./LD_stack/LD.clump.result.eb_test_GA_",i1,".rdata"))
  LD.result.eb.stack.test = LD.result.list[[1]]
  
  
  LD.result.stack.com <- rbind(LD.result.stack,LD.result.2D.stack,
                               LD.result.eb.stack,
                               LD.result.eb.stack.test) %>% 
    filter(eth.vec!="EUR")
  
  sample_size =  as.character(LD.result.stack.com$method_vec)
  ssp = as.character(c("15000","45000","80000","100000"))
  cau_vec <- as.character(LD.result.stack.com$l_vec)
  csp <- c(0.01,0.001,0.0005)
  for(l in 1:3){
    idx <- which(LD.result.stack.com$l_vec==l)
    cau_vec[idx] <- paste0("Causal SNPs Proportion = ",csp[l])
  }
  cau_vec = factor(cau_vec,levels = paste0("Causal SNPs Proportion = ",csp))
  for(m in 1:1){
    idx <- which(LD.result.stack.com$m_vec==m)
    sample_size[idx] <- ssp[m]
  }
  sample_size = factor(sample_size,
                       levels = c("15000","45000","80000","100000"))
  LD.result.stack.com <- cbind(LD.result.stack.com,cau_vec,sample_size)  
    LD.result.stack.sub = LD.result.stack.com %>% 
    select(colnames(LD.clump.result.sub))
    LD.clump.result.plot.sub <- rbind(LD.clump.result.sub,LD.result.stack.sub) %>% 
      filter(method_vec!="P+T") 
    
    jdx <- which(LD.clump.result.plot.sub$method_vec=="Best EUR PRS")
    LD.clump.result.plot.sub$method_vec[jdx] = "Best EUR PRS (C+T)"
    jdx <- which(LD.clump.result.plot.sub$method_vec=="Best EUR SNP + target coefficients")
    LD.clump.result.plot.sub$method_vec[jdx] = "Best EUR SNP + target coefficients (C+T)"
    jdx <- which(LD.clump.result.plot.sub$method_vec=="Best EUR SNP + EB")
    LD.clump.result.plot.sub$method_vec[jdx] = "Best EUR SNP + EB(C+T)"
    method_vec = factor(LD.clump.result.plot.sub$method_vec,
                        levels = c("C+T","C+T max","C+T SL",
                                   "LDpred2",
                                   
                                   "Best EUR PRS (C+T)",
                                   "Best EUR SNP + target coefficients (C+T)",
                                   "Best EUR SNP + EB(C+T)",
                                   "EUR LDpred2",
                                   "2DLD",
                                   "2DLD-max",
                                   "2DLD-SL",
                                   "2DLD-eb",
                                   "2DLD-max-eb",
                                   "2DLD-SL-eb",
                                   "2DLD-SL-eb-test",
                                   "MEBayes"))
  
    LD.clump.result.plot.sub$method_vec = method_vec
    
    
    LD.clump.result.plot.sub.1 = LD.clump.result.plot.sub %>% 
      filter(method_vec %in%c("C+T","C+T max","C+T SL",
                          "LDpred2",
                          
                          "Best EUR PRS (C+T)",
                          "Best EUR SNP + target coefficients (C+T)",
                          "Best EUR SNP + EB(C+T)",
                          "EUR LDpred2"))
    
    
    
    p <- ggplot(LD.clump.result.plot.sub.1,aes(x= sample_size,y=r2.vec,group=method_vec))+
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
      scale_fill_manual(values = getPalette(colourCount)) +
      theme(axis.text = element_text(size = rel(0.9)),
            legend.text = element_text(size = rel(0.9)))+
      ggtitle("Prediction performance comparasion (Sample Size for EUR = 100k)")
    p
    png(file = paste0("./LD_stack/sub1_method_compare_result_size_",m,"_summary_GA_",i1,".png"),
        width = 13, height = 8, res = 300,units = "in")
    print(p)
    dev.off()
    
    LD.clump.result.plot.sub.2 = LD.clump.result.plot.sub %>% 
      filter(method_vec %in%c("Best EUR SNP + EB(C+T)",
                              "EUR LDpred2",
                              "2DLD",
                              "2DLD-max",
                              "2DLD-SL",
                              "2DLD-eb",
                              "2DLD-max-eb",
                              "2DLD-SL-eb",
                              "2DLD-SL-eb-test",
                              "MEBayes"))
    
    p <- ggplot(LD.clump.result.plot.sub.2,aes(x= sample_size,y=r2.vec,group=method_vec))+
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
      scale_fill_manual(values = getPalette(colourCount)) +
      theme(axis.text = element_text(size = rel(0.9)),
            legend.text = element_text(size = rel(0.9)))+
      ggtitle("Prediction performance comparasion (Sample Size for EUR = 100k)")
    p
    LD.clump.result.plot.sub.2 %>% 
      filter(method_vec%in%c("2DLD-SL-eb-test"))
    png(file = paste0("./LD_stack/sub2_method_compare_result_size_",m,"_summary_GA_",i1,".png"),
        width = 13, height = 8, res = 300,units = "in")
    print(p)
    dev.off()
    
    
    
    
    
  }
  LD.clump.result.15k <- LD.clump.result %>% filter(m_vec==1) %>% 
    select(eth.vec,r2.vec,cau_vec,sample_size,method_vec) %>% 
    mutate(method_vec = as.character(method_vec))
  method_vec = as.character(LD.clump.result.15k$method_vec)
  method_vec = factor(method_vec,
                      
                      levels = c("P+T",
                                 "Best EUR PRS",
                                 "Best EUR SNP + target coefficients",
                                 "Best EUR SNP + EB",
                                 "2DLD",
                                 "2WLD"))
  LD.clump.result.15k$method_vec = method_vec
  
  load("R2.simulation.hapmap3.RData")
  res$cau_vec = as.character(res$cau_vec)
  idx <- which(res$cau_vec=="Causal SNPs Proportion = 5e-4")
  res$cau_vec[idx] = "Causal SNPs Proportion = 5e-04"
  LD.result.LDpred <- res %>% 
    rename(eth.vec = race,
           r2.vec = R2,
           method_vec = Method) %>% 
    filter(GA==i1) %>% 
    select(eth.vec,r2.vec,cau_vec,sample_size,method_vec)
  LD.clump.result.15k <- rbind(LD.clump.result.15k,LD.result.LDpred)
  method_vec = as.character(LD.clump.result.15k$method_vec)
  method_vec = factor(method_vec,
                      levels = c("P+T","LDpred2",
                                 "Best EUR PRS",
                                 "Best EUR SNP + target coefficients",
                                 "Best EUR SNP + EB",
                                 "EUR LDpred2",
                                 "2DLD",
                                 "2WLD",
                                 "MEBayes"))
  LD.clump.result.15k$method_vec = method_vec
  colourCount = 9
  getPalette = colorRampPalette(brewer.pal(9, "Paired"))
  
  p <- ggplot(LD.clump.result.15k,aes(x= sample_size,y=r2.vec,group=method_vec))+
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
    scale_fill_manual(values = getPalette(colourCount)) +
    theme(axis.text = element_text(size = rel(0.9)),
          legend.text = element_text(size = rel(0.9)))+
    ggtitle("Prediction performance comparasion (Sample Size for EUR = 100k)")
  p
  png(file = paste0("./all_method_compare_result_size_",1,"_summary_GA_",i1,".png"),
      width = 13, height = 8, res = 300,units = "in")
  print(p)
  dev.off()
  
}




setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA")
library(ggplot2)
library(ggsci)
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
p
png(file = paste0("./heritability.png"),
    width = 10, height = 8, res = 300,units = "in")
p
dev.off()





#load LD clump results
load("LD.clump.result.rdata")
LD.clump.result <- LD.result.list[[1]]

sample_size <- factor(rep(c("15000","45000","80000","100000"),4*3*2),
                      levels=c("15000","45000","80000","100000"))

cau_vec <- as.character(LD.clump.result$l_vec)
csp <- c(0.01,0.001,0.0005)
LD.clump.result <- LD.clump.result[LD.clump.result$eth.vec!="EUR",]
method_vec <- rep("P+T",nrow(LD.clump.result))
LD.clump.result$method_vec = method_vec
LD.clump.result.PT <- LD.clump.result
#Best EUR result
load("best_eur_snp_result.rdata")
method <- c("eurcoef","tarcoef","eb")
method_nameupdate <- c("Best EUR PRS","Best EUR SNP + target coefficients","Best EUR SNP + EB")
for(q in 1:length(method)){
  idx <- which(LD.clump.result$method_vec==method[q])
  LD.clump.result$method_vec[idx] <-   method_nameupdate[q]
}
LD.clump.result.EUR <- LD.clump.result
# 
# 

# 
load("LD.clump.result.two.dim.rdata")

LD.clump.result <- LD.result.list[[1]]

sample_size <- factor(rep(c("15000","45000","80000","100000"),4*3*2),
                      levels=c("15000","45000","80000","100000"))

cau_vec <- as.character(LD.clump.result$l_vec)
csp <- c(0.01,0.001,0.0005)
#LD.clump.result <- LD.clump.result[LD.clump.result$eth.vec!="EUR",]
method_vec <- rep("2DLD",nrow(LD.clump.result))
LD.clump.result$method_vec = method_vec
LD.clump.result.L2 <- LD.clump.result



load("LD.clump.result.eb.rdata")

LD.clump.result <- LD.result.list[[1]]

sample_size <- factor(rep(c("15000","45000","80000","100000"),4*3*2),
                      levels=c("15000","45000","80000","100000"))

cau_vec <- as.character(LD.clump.result$l_vec)
csp <- c(0.01,0.001,0.0005)
#LD.clump.result <- LD.clump.result[LD.clump.result$eth.vec!="EUR",]
method_vec <- rep("2DLD-EB",nrow(LD.clump.result))
LD.clump.result$method_vec = method_vec
LD.clump.result.L2EB <- LD.clump.result

LD.clump.result <- rbind(LD.clump.result.L2EB,LD.clump.result.L2,LD.clump.result.EUR,LD.clump.result.PT)





#LD.clump.result <- rbind(LD.clump.result.L2,LD.clump.result.EUR,LD.clump.result.PT)




for(l in 1:3){
  idx <- which(LD.clump.result$l_vec==l)
  cau_vec[idx] <- paste0("Causal SNPs Proportion = ",csp[l])
}
cau_vec <- factor(cau_vec,levels = paste0("Causal SNPs Proportion = ",csp))
sample_size <- as.character(LD.clump.result$m_vec)
sample_option = c("15000","45000","80000","100000")
for(m in 1:4){
  idx <- which(LD.clump.result$m_vec==m)
  sample_size[idx] <- paste0(sample_option[m])
}

sample_size = factor(sample_size,levels=c("15000","45000","80000","100000"))

LD.clump.result <- cbind(LD.clump.result,cau_vec,sample_size)

LD.clump.result$method_vec <- factor(LD.clump.result$method_vec,
                                     levels = c("P+T",
                                                "Best EUR PRS",
                                                "Best EUR SNP + target coefficients",
                                                "Best EUR SNP + EB",
                                                "2DLD","2DLD-EB"))


p <- ggplot(LD.clump.result,aes(x= sample_size,y=r2.vec,group=method_vec))+
  geom_bar(aes(fill=method_vec),
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
  ggtitle("Prediction performance comparasion (Sample Size for EUR = 100k)")
p
png(file = paste0("./method_compare_result_summary.png"),
    width = 13, height = 8, res = 300,units = "in")
p
dev.off()

idx <- which(LD.clump.result$m_vec==1)
p <- ggplot(LD.clump.result[idx,],aes(x= sample_size,y=r2.vec,group=method_vec))+
  geom_bar(aes(fill=method_vec),
           stat="identity",
           position = position_dodge())+
  #geom_point(aes(color=method_vec))+
  theme_Publication()+
  ylab("R2")+
  xlab("Sample Size")+
  labs(fill = "Method")+
  facet_grid(vars(cau_vec),vars(eth.vec))+
  scale_fill_nejm()+
  ggtitle("Prediction comparasion (Sample Size for EUR = 100k)")
p



png(file = paste0("./method_compare_result_summary_15000.png"),
    width = 12, height = 8, res = 300,units = "in")
p
dev.off()



#load Jin's results
load("results.hapmap3.RData")
res$cau_vec <- as.character(res$cau_vec)
idx <- which(res$cau_vec=="Causal SNPs Proportion = 5e-4")
res$cau_vec[idx] <- "Causal SNPs Proportion = 5e-04"
res <- res[,-6]
colnames(res) <- colnames(LD.clump.result)
res$method_vec <- paste0(res$method_vec," (hap3)")
LD.clump.result <- rbind(LD.clump.result,res)
LD.clump.result$method_vec <- factor(LD.clump.result$method_vec,
                                     levels = c("P+T",
                                                "P+T (hap3)",
                                                "Best EUR PRS",
                                                "EUR P+T (hap3)",
                                                "Best EUR SNP + target coefficients",
                                                "Best EUR SNP + EB",
                                                "LDpred (hap3)",
                                                "EUR LDpred (hap3)",
                                                "ME-Bayes (hap3)",
                                                "2DLD",
                                                "2DLD-EB"))
idx <- which(LD.clump.result$m_vec==1&
               LD.clump.result$eth.vec!="SAS")
library(RColorBrewer)
colourCount = length(unique(LD.clump.result$method_vec))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
p <- ggplot(LD.clump.result[idx,],aes(x= sample_size,y=r2.vec,group=method_vec))+
  geom_bar(aes(fill=method_vec),
           stat="identity",
           position = position_dodge())+
  #geom_point(aes(color=method_vec))+
  theme_Publication()+
  ylab("R2")+
  xlab("Sample Size")+
  labs(fill = "Method")+
  facet_grid(vars(cau_vec),vars(eth.vec))+
  scale_fill_manual(values = getPalette(colourCount))+
  ggtitle("Prediction comparasion (Sample Size for EUR = 100k)")
p
png(file = paste0("./all_method_compare_result_summary_15000.png"),
    width = 13, height = 8, res = 300,units = "in")
p
dev.off()


idx <- which(LD.clump.result$m_vec==2&
               LD.clump.result$eth.vec!="SAS")
p <- ggplot(LD.clump.result[idx,],aes(x= sample_size,y=r2.vec,group=method_vec))+
  geom_bar(aes(fill=method_vec),
           stat="identity",
           position = position_dodge())+
  #geom_point(aes(color=method_vec))+
  theme_Publication()+
  ylab("R2")+
  xlab("Sample Size")+
  labs(fill = "Method")+
  facet_grid(vars(cau_vec),vars(eth.vec))+
  scale_fill_manual(values = getPalette(colourCount))+
  ggtitle("Prediction comparasion (Sample Size for EUR = 100k)")
p
png(file = paste0("./all_method_compare_result_summary_15000.png"),
    width = 13, height = 8, res = 300,units = "in")
p
dev.off()
