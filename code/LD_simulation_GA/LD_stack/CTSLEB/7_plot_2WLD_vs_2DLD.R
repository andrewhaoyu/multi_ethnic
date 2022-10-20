setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA")
source("../../code/LD_simulation_large/theme_Publication.R")
library(ggplot2)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(grid)
library(gridExtra)

#
#2DLD
load(paste0("LD.clump.result.two.dim_GA_",i1,".rdata"))
LD.clump.result <- LD.result.list[[1]]
method_vec <- rep("2DLD",nrow(LD.clump.result))
LD.clump.result$method_vec = method_vec
LD.clump.result.2DLD = LD.clump.result

#2DLD + EB
# load(paste0("LD.clump.result.eb_GA_",i1,".rdata"))
# LD.clump.result <- LD.result.list[[1]]
# method_vec <- rep("2DLD + EB coefficients",nrow(LD.clump.result))
# LD.clump.result$method_vec = method_vec
# LD.clump.result.EB = LD.clump.result
#2WLD
load(paste0("LD.clump.result.two_way_GA_",i1,".rdata"))
LD.clump.result <- LD.result.list[[1]]
method_vec <- rep("2WLD",nrow(LD.clump.result))
LD.clump.result$method_vec = method_vec
LD.clump.result.2WLD = LD.clump.result

LD.clump.result <- rbind(LD.clump.result.2DLD,LD.clump.result.2WLD)
LD.clump.result = LD.clump.result %>% 
  mutate(method_vec_update = case_when(
    method_vec == "2DLD" ~ "Alternative way",
    method_vec == "2WLD" ~ "Proposed way"
  ))
LD.clump.result$method_vec_update <- factor(LD.clump.result$method_vec_update,
                                     levels = c("Alternative way",
                                                "Proposed way"))
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
for(m in 1:4){
  LD.clump.result.sub <- LD.clump.result %>% filter(m_vec==m) 
 
  p <- ggplot(LD.clump.result.sub,aes(x= sample_size,y=r2.vec,group=method_vec_update))+
    geom_bar(aes(fill=method_vec_update),
             stat="identity",
             position = position_dodge())+
    #geom_point(aes(color=method_vec))+
    theme_Publication()+
    ylab("R2")+
    xlab("Sample Size")+
    labs(fill = "Method")+
    facet_grid(vars(cau_vec),vars(eth.vec))+
    #scale_fill_nejm()+
    scale_fill_Publication() +
    theme(axis.text = element_text(size = rel(0.9)),
          legend.text = element_text(size = rel(0.9)))+
    ggtitle("Fixed common SNPs heritability with strong negative selection")
  p
  png(file = paste0("./LD_stack/simulation_2wld_2dld/method_compare_result_size_",m,"_summary_GA_",i1,".png"),
      width = 19, height = 12, res = 300,units = "in")
  print(p)
  dev.off()
  
}
