setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_new")
library(ggplot2)
library(ggsci)
#standard LD clumping results
load("LD.clump.result.rdata")
LD.clump.result <- LD.result.list[[1]]
sample_size <- factor(rep(c("15000","45000","80000","100000"),15),
                      levels=c("15000","45000","80000","100000"))

cau_vec <- as.character(LD.clump.result$l_vec)
csp <- c(0.01,0.001,0.0005)
for(l in 1:3){
  idx <- which(LD.clump.result$l_vec==l)
  cau_vec[idx] <- paste0("Causal SNPs Proportion = ",csp[l])
}
cau_vec <- factor(cau_vec,
                  levels = paste0("Causal SNPs Proportion = ",csp))

LD.clump.result <- cbind(LD.clump.result,sample_size,cau_vec)
save(LD.clump.result,file = "LD.clump.result_090420_P+T.rdata")
load("LD.clump.result_090420_P+T.rdata")


LD.clump.result$eth.vec <- factor(LD.clump.result$eth.vec,
                                  levels =c("EUR","AFR","AMR","EAS","SAS"))
p <- ggplot(LD.clump.result,aes(x= sample_size,y=r2.vec,group=eth.vec))+
  geom_line(aes(color=eth.vec))+
  geom_point(aes(color=eth.vec))+
  theme_Publication()+
  ylab("R2")+
  xlab("Training sample size")+
  labs(color = "Ethnic group")+
  facet_grid(cols=vars(cau_vec))+
  scale_color_nejm()+
  ggtitle("Prediction performance across ethnic groups using P+T")
p
png(file = paste0("./LD_clumping_result_summary.png"),
    width = 10, height = 8, res = 300,units = "in")
p
dev.off()

#genearte plot by varying the p-value threshold
LD.clump.result <- LD.result.list[[2]]
sample_size  <- LD.clump.result$m_vec
sample_size_option <- paste0( "Sample Size = ", c("15000","45000","80000","100000"))
for(m in 1:4){
  idx <- which(LD.clump.result$m_vec==m)
  sample_size[idx] <- sample_size_option[m]
}
sample_size <- factor(sample_size,
                      levels=paste0("Sample Size = ", c("15000","45000","80000","100000")))

cau_vec <- as.character(LD.clump.result$l_vec)
csp <- c(0.01,0.001,0.0005)
for(l in 1:3){
  idx <- which(LD.clump.result$l_vec==l)
  cau_vec[idx] <- paste0("Causal SNPs Proportion = ",csp[l])
}
cau_vec <- factor(cau_vec,
                  levels = paste0("Causal SNPs Proportion = ",csp))
LD.clump.result <- cbind(LD.clump.result,sample_size,cau_vec)



  
  p <- ggplot(LD.clump.result,aes(x= log10(pthres.vec),y=r2.vec,group=eth.vec))+
    geom_line(aes(color=eth.vec))+
    geom_point(aes(color=eth.vec))+
    theme_Publication()+
    ylab("R2")+
    xlab("log10(P-value)")+
    labs(color = "Ethnic group")+
    facet_grid(vars(cau_vec),vars(sample_size))+
    scale_color_nejm()
  p
  png(file = paste0("./LD_clumping_result_p_thres.png"),
      width = 15, height = 10, res = 300,units = "in")
  p
  dev.off()

  
#Best EUR result
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
  LD.clump.result.plot <- LD.clump.result[idx,]
  jdx <- which(LD.clump.result.plot$method_vec%in%
                 c("P+T","Best EUR PRS","Best EUR SNP + target coefficients",
                   "Best EUR SNP + EB","LDpred (hap3)","EUR LDpred (hap3)","ME-Bayes (hap3)","2DLD","2DLD-EB"))
  LD.clump.result.plot <- LD.clump.result.plot[jdx,]
  
  #update some of the names
  method <- as.character(LD.clump.result.plot$method_vec)
  kdx <- which(method =="LDpred (hap3)")
  method[kdx] <- "LDpred"
  kdx <- which(method =="EUR LDpred (hap3)")
  method[kdx] <- "EUR LDpred PRS"
  kdx <- which(method =="ME-Bayes (hap3)")
  method[kdx] <- "ME-Bayes"
  kdx <- which(method =="Best EUR PRS")
  method[kdx] <- "EUR P+T PRS"
  kdx <- which(method =="Best EUR SNP + target coefficients")
  method[kdx] <- "SNPs in EUR P+T + target coefficients"
  kdx <- which(method =="Best EUR SNP + EB")
  method[kdx] <- "SNPs in EUR P+T + EB coefficients"
  method.fac <- factor(method,
                       levels = c("P+T",
                                  "LDpred",
                                  "EUR P+T PRS",
                                  "EUR LDpred PRS",
                                  "SNPs in EUR P+T + target coefficients",
                                  "SNPs in EUR P+T + EB coefficients",
                                  "ME-Bayes",
                                  "2DLD",
                                  "2DLD-EB"))
  LD.clump.result.plot$method_vec <- method.fac
  
  
  save(LD.clump.result.plot,file ="data_for_plot.rdata")
  p <- ggplot(LD.clump.result.plot,aes(x= sample_size,y=r2.vec,group=method_vec))+
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
  