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


p <- ggplot(LD.clump.result,aes(x= sample_size,y=r2.vec,group=eth.vec))+
  geom_line(aes(color=eth.vec))+
  geom_point(aes(color=eth.vec))+
  theme_Publication()+
  ylab("R2")+
  xlab("Training sample size")+
  guides(color=guide_legend(title="Ethnic group"))+
  facet_grid(cols=vars(cau_vec))+
  scale_color_nejm()
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
    guides(color=guide_legend(title="Ethnic group"))+
    facet_grid(vars(cau_vec),vars(sample_size))+
    scale_color_nejm()
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
  
  LD.clump.result <- rbind(LD.clump.result,LD.clump.result.PT)
  for(l in 1:3){
    idx <- which(LD.clump.result$l_vec==l)
    cau_vec[idx] <- paste0("Causal SNPs Proportion = ",csp[l])
  }
  cau_vec <- factor(cau_vec,
                    levels = paste0("Causal SNPs Proportion = ",csp))
  
  LD.clump.result <- cbind(LD.clump.result,sample_size,cau_vec)
  
  
  p <- ggplot(LD.clump.result,aes(x= sample_size,y=r2.vec,group=method_vec))+
    geom_line(aes(color=method_vec))+
    geom_point(aes(color=method_vec))+
    theme_Publication()+
    ylab("R2")+
    xlab("Sample Size")+
    guides(color=guide_legend(title="Method"))+
    facet_grid(vars(cau_vec),vars(eth.vec))+
    scale_color_nejm()
  png(file = paste0("./method_compare_result_summary.png"),
      width = 12, height = 8, res = 300,units = "in")
  p
  dev.off()
  
  
  