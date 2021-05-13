setwd("/Users/zhangh24/GoogleDrive/multi_ethnic/result/23andme/PT_summary/summary")
eth_group = c("european","african_american",
              "latino","east_asian","south_asian")
eth_name = c("European","African American",
             "Latino","East Asian","South Asian")
eth <- c("EUR","AFR","AMR","EAS","SAS")
bin_trait <- c("any_cvd","depression",
               "iqb.sing_back_musical_note",
               "migraine_diagnosis",
               "morning_person")
con_trait <- c(       "heart_metabolic_disease_burden",
                      "height"
)
trait = c("any_cvd","depression",
          "heart_metabolic_disease_burden",
          "height",
          "iqb.sing_back_musical_note",
          "migraine_diagnosis",
          "morning_person")
trait_name = c("Any CVD","Depression",
               "Heart metabolic disease burden",
               "Height",
               "SBMN",
               "Migraine Diagnosis",
               "Morning Person")
method_vec = c("PT","tdld")
method_name = c("P + T","TDLD")
source("/Users/zhangh24/GoogleDrive/multi_ethnic/code/stratch/theme_publication.R")
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
method = "PT"
#load all results
plot.data.list = list()
library(data.table)
temp = 1
for(i1 in 1:length(method_vec)){
  for(l in 1:7){
    for(i in 2:5){
      result = read.csv(paste0(eth_group[i],"_",trait[l],"_",method_vec[i1]))
      
      plot.data = data.frame(result = as.numeric(max(result)),
                             eth = eth_name[i],
                             trait = trait_name[l],
                             Method = method_name[i1])
      
      plot.data.list[[temp]] = plot.data
      temp = temp+1
    }
  }
}
plot.data = rbindlist(plot.data.list)
plot.data$index = rep("1",nrow(plot.data))
plot.data$eth = factor(plot.data$eth,
                       levels = c("European","African American",
                                  "Latino","East Asian","South Asian"))
plot.data.sub = plot.data %>% 
  filter(trait %in% c("Any CVD","Depression",
                      "SBMN",
                      "Migraine Diagnosis",
                      "Morning Person")) 

p = ggplot(plot.data.sub)+
  geom_bar(aes(x = index,y = result,fill=Method),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  ylab("AUC")+
  facet_grid(vars(trait),vars(eth))+
  theme_Publication()+
  scale_fill_Publication()+
  coord_cartesian(ylim = c(0.5, 0.64)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#library(RColorBrewer)

png(filename = "../../comparasion_plot/bin_comparasion.png",width=12,height = 8,units = "in",res = 300)
print(p)
dev.off()

plot.data.sub = plot.data %>% 
  filter(trait %in% c("Heart metabolic disease burden","Height")) 


p = ggplot(plot.data.sub)+
  geom_bar(aes(x = index,y = result,fill=Method),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  ylab("R2")+
  facet_grid(vars(trait),vars(eth),scales = "free")+
  #facet_wrap(vars(eth))+
  theme_Publication()+
  scale_fill_Publication()+
  #coord_cartesian(ylim = c(0.5, 0.64)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  #theme(strip.text.y = element_blank())
#library(RColorBrewer)
png(filename = "../../comparasion_plot/con_comparasion.png",width=10,height = 8,units = "in",res = 300)
print(p)
dev.off()















load(paste0("LD.clump.result.CT.rdata"))
load(paste0("LD.clump.result.SCT.rdata"))
load(paste0("eur.snp.reult.rdata"))
load(paste0("weightedprs.result.rdata"))
load(paste0("LD.clump.result.2DLD.rdata"))
load(paste0("LD.clump.result.EB.rdata"))
load(paste0("LD.clump.result.alleth.EB.rdata"))
load(paste0("LDpred2.result.rdata"))
load(paste0("LDpredEUR.result.rdata"))

LD.clump.result <- LD.result.list[[1]] %>% 
  mutate(method_vec = rep("C+T"))

TDLD.result = TDLD.result %>% 
  filter(method_vec=="TDLD")

alleth.EB.result = alleth.EB.result %>% 
  filter(method_vec=="TDLD-SLEB (all ethnics)")
weightedprs.result = weightedprs.result %>% 
  mutate(method_vec = "Weighted PRS")
prediction.result <- rbind(LD.clump.result,
                           SCT.clump.result,
                           LDpred2.result,
                           eursnp.result,
                           LDpredEUR.result,
                           weightedprs.result,
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
                                        "TDLD",
                                        "TDLD-EB",
                                        "TDLD-SLEB",
                                        "TDLD-SLEB (all ethnics)"
                             ))) %>% 
  mutate(ga_arc = case_when(ga_vec==1 ~"Fixed common SNP heritability",
                            ga_vec==2 ~"Strong negative selection",
                            ga_vec==3 ~"Less correlation",
                            ga_vec==4 ~"No negative selection",
                            ga_vec==5 ~"Mild negative selection",
  ))
uvals = unique(prediction.result$method_vec)

n.single = 9


single.color =  brewer.pal(n.single, "Blues")[c(4,5,7)]
n.EUR = 9


EUR.color = brewer.pal(n.EUR, "Greens")[c(4,5,6,7)]


n.multi = 9
multi.color = brewer.pal(n.multi, "Oranges")[c(3:7)]
colour = c(single.color,EUR.color,multi.color)
col_df = tibble(
  colour = c(single.color,EUR.color,multi.color),
  method_vec = uvals,
  category = case_when(method_vec%in%c("C+T",
                                       "SCT",
                                       "LDpred2") ~ "Single ethnic method",
                       method_vec%in%c("Best EUR SNP (C+T)",
                                       "Best EUR SNP + target coefficients (C+T)",
                                       "Best EUR SNP + EB coefficients (C+T)",
                                       "Best EUR PRS (LDpred2)"
                       ) ~ "EUR PRS based method",
                       method_vec%in%c("Weighted PRS",
                                       "TDLD",
                                       "TDLD-EB",
                                       "TDLD-SLEB",
                                       "TDLD-SLEB (all ethnics)") ~ "Multi ethnic method")
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
m = 1
for(m in 1:4){
  for(i1 in 1:5){
    
    prediction.result.sub <- prediction.result %>% 
      filter(ga_vec==i1&
               eth.vec!="EUR"&
               m_vec ==m)
    title = prediction.result.sub$ga_arc[1]
    legs = lapply(sort(unique(col_df$category)), run_plot)
    
    legs = lapply(legs, getLegend)
    p.leg = plot_grid(NULL,NULL,legs[[1]],legs[[2]], legs[[3]],NULL,align="v",ncol=1,rel_heights=c(1,1,0.8,0.9,1.2,1.2))
    print(p.leg)
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
      ggtitle(title)+
      theme(legend.position = "none")
    p = plot_grid(p.null,p.leg,nrow=1,rel_widths = c(3.5,1))
    print(p)
    png(file = paste0("./method_compare_result_size_",m,"_summary_GA_",i1,".png"),
        width = 13, height = 8, res = 300,units = "in")
    print(p)
    dev.off()
    
    
  }
  
  
}






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
