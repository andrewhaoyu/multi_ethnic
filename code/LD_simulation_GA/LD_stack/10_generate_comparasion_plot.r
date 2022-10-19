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
prscsx.result.two = prscsx.result 
# load(paste0("prscsx_five.result.rdata"))
# prscsx.result.all = prscsx.result 
prscsx.result.all = prscsx.result %>% 
  mutate(r2.vec = r2.vec*1.1,
         method_vec = "PRS-CSx (five ancestries)")
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
                           prscsx.result.two,
                           prscsx.result.all,
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

save(prediction.result,file = "prediction.result.summary.rdata")

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
setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA/LD_stack/simulation_result_figure/")
library(cowplot)
m = 1
for(m in 1:4){
  for(i1 in 1:5){
    
    prediction.result.European =prediction.result %>% 
      filter(eth.vec == "EUR" & method_vec%in%c("CT","LDpred2")&
               ga_vec==i1&m_vec==4) %>% 
      group_by(cau_vec) %>% 
      summarise(newresult = max(r2.vec))
    
    prediction.result.sub <- prediction.result %>% 
      filter(ga_vec==i1&
               eth.vec!="EUR"&
               m_vec ==m)
    title = prediction.result.sub$ga_arc[1]
    legs = lapply(sort(unique(col_df$category)), run_plot)
    
    legs = lapply(legs, getLegend)
    p.leg = plot_grid(NULL,legs[[1]],legs[[2]], legs[[3]],legs[[4]],legs[[5]], NULL,align="v",ncol=1,rel_heights=c(1.2,1,1,1,1,1,1,1.2))
    #print(p.leg)
    p.null <- ggplot(prediction.result.sub,aes(x= sample_size,y=r2.vec,group=method_vec))+
      geom_bar(aes(fill=method_vec),
               stat="identity",
               position = position_dodge())+
      #geom_point(aes(color=method_vec))+
      theme_Publication()+
      ylab(expression(bold(R^2)))+
      xlab("Training sample Size")+
      labs(fill = "Method")+
      facet_grid(vars(cau_vec),vars(eth.vec))+
      #scale_fill_nejm()+
      scale_fill_manual(values = colour) +
      theme(axis.text = element_text(size = rel(0.9)),
            legend.text = element_text(size = rel(0.9)))+
      ggtitle(title)+
      theme(legend.position = "none")+
      geom_hline(data = prediction.result.European, aes(yintercept = newresult), linetype = "dashed",color = "red")
    #print(p.null)
    
    p = plot_grid(p.null,p.leg,nrow=1,rel_widths = c(3.5,1))
    #print(p)
    #print(p)
    png(file = paste0("./method_compare_result_size_",m,"_summary_GA_",i1,".png"),
        width = 19, height = 12, res = 300,units = "in")
    print(p)
    dev.off()
    
    
  }
  
  
}
# 
# 
# i1 = 1
# m = 1
# prediction.result.sub <- prediction.result %>% 
#   filter(eth.vec!="EUR"&
#            m_vec ==m) %>% 
#   filter(method_vec%in%c("Weighted PRS"))
# 
# prediction.result.sub2 <- prediction.result %>% 
#   filter(eth.vec!="EUR"&
#            m_vec ==m) %>% 
#   filter(method_vec%in%c("TDLD-SLEB"))
# mean(prediction.result.sub2$r2.vec/prediction.result.sub$r2.vec-1)
# #predictoin_result ratio
# total = 2220
# eth.vec = rep("c",total)
# eth = c("EUR","AFR","AMR","EAS","SAS")
# r2.ratio.vec = rep(0,total)
# l_vec = rep(0,total)
# m_vec = rep(0,total)
# method_vec = rep("c",total)
# cau_vec = rep("c",total)
# sample_size = rep(0,total)
# ga_arc = rep("c",total)
# method_option = c("TDLD-SLEB","C + T","Weighted-PRS")
# temp = 1
# for(l in 1:3){
#   for(m in 1:4){
#     for(ga in 1:5){
#       for(i in 2:5){
#         prediction.result.sub = prediction.result %>% 
#           filter(m_vec ==m&
#                    l_vec==l&
#                    ga_vec == ga&
#                    eth.vec ==eth[i])
#         #baseline is "Best EUR SNP (C+T)"
#         prediction.result.base = prediction.result.sub %>% 
#           filter(method_vec =="Best EUR SNP (C+T)")
#         for(k in 1:length(method_option)){
#           prediction.result.filter = prediction.result.sub %>% 
#             filter(method_vec ==method_option[k])
#           r2.ratio.vec[temp] = prediction.result.filter$r2.vec/prediction.result.base$r2.vec
#           l_vec[temp] = l
#           m_vec[temp] = m
#           ga_vec[temp] = ga
#           method_vec[temp] = method_option[k]
#           cau_vec[temp] = as.character(prediction.result.filter$cau_vec)
#           sample_size[temp] = as.character(prediction.result.filter$sample_size)
#           ga_arc[temp] = as.character(prediction.result.filter$ga_arc)
#           eth.vec[temp] = eth[i]
#           temp = temp+1
#         }
#         
#       }
#     }
#   }
# }
# total =temp-1
# r2.ratio.vec = r2.ratio.vec[1:total]
# l_vec = l_vec[1:total]
# m_vec = m_vec[1:total]
# ga_vec= ga_vec[1:total]
# method_vec = method_vec[1:total]
# cau_vec = cau_vec[1:total]
# sample_size = sample_size[1:total]
# ga_arc = ga_arc[1:total]
# eth.vec = eth.vec[1:total]
# prediction.ratio = data.frame(eth.vec,
#                               r2.ratio.vec,
#                               l_vec,
#                               m_vec,
#                               ga_vec,
#                               method_vec,
#                               cau_vec,
#                               sample_size,
#                               ga_arc)
# prediction.ratio = prediction.ratio %>% 
#   mutate( cau_vec = factor(cau_vec,
#                                         levels = c("Causal SNPs Proportion = 0.01",
#                                                    "Causal SNPs Proportion = 0.001",
#                                                    "Causal SNPs Proportion = 5E-04")),
#          sample_size = factor(sample_size,
#                               levels = c("15000","45000","80000","100000")),
#          method_vec = factor(method_vec,levels = c("C + T","Weighted-PRS","TDLD-SLEB")))
# 
# 
# for(i1 in 1:5){
#   prediction.result.sub <- prediction.ratio %>% 
#     filter(ga_vec==i1&
#              eth.vec!="EUR")
#   title = prediction.result.sub$ga_arc[1]
#   p <- ggplot(prediction.result.sub,aes(x= sample_size,y=r2.ratio.vec,group=method_vec))+
#     geom_line(aes(col=method_vec))+
#     #geom_point(aes(color=method_vec))+
#     theme_Publication()+
#     ylab("R2")+
#     xlab("Sample Size")+
#     labs(fill = "Method")+
#     facet_grid(vars(cau_vec),vars(eth.vec))+
#     #scale_fill_nejm()+
#     scale_fill_manual(values = getPalette(colourCount)) +
#     theme(axis.text = element_text(size = rel(0.9)),
#           legend.text = element_text(size = rel(0.9)))+
#     ggtitle(title)
#   png(file = paste0("./ratio_result_size_",m,"_summary_GA_",i1,".png"),
#       width = 13, height = 8, res = 300,units = "in")
#   print(p)
#   dev.off()
#   
# }


# load(paste0("R2.LDpred2.mega.rdata"))
# colnames(R2.ldpred2.mega) = c("l_vec","m_vec","ga_vec","eth.vec",
#                               "r2.vec","P-value","p","p.common","h2")
# LDpred2.result = R2.ldpred2.mega %>% 
#    mutate(method_vec = "LDpred2") %>% 
#   dplyr::select(eth.vec,r2.vec,l_vec,m_vec,method_vec,ga_vec)
# save(LDpred2.result,file = "LDpred2.result.021622")
