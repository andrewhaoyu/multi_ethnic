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
               "Sing back musical note",
               "Migraine Diagnosis",
               "Morning Person")
source("/Users/zhangh24/GoogleDrive/multi_ethnic/code/stratch/theme_publication.R")
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
method = "PT"
#load all results
plot.data.list = list()
library(data.table)
temp = 1
  for(l in 1:7){
    for(i in 1:5){
    result = read.csv(paste0(eth_group[i],"_",trait[l],"_",method))
    
    plot.data = data.frame(pthres = pthres[(1+length(pthres)-length(result)):length(pthres)],
                           result = as.numeric((result)),
                           eth = eth_name[i],
                           trait = trait_name[l])
    
    plot.data.list[[temp]] = plot.data
    temp = temp+1
  }
  }
plot.data = rbindlist(plot.data.list)
plot.data$eth = factor(plot.data$eth,
                       levels = c("European","African American",
                                 "Latino","East Asian","South Asian"))

plot.data.sub = plot.data %>% 
  filter(trait %in% c("Any CVD","Depression",
                      "Sing back musical note",
                      "Migraine Diagnosis",
                      "Morning Person"))

# result = read.csv(paste0(eth_group[i],"_",trait[l],"_",method))
# library(tidyverse)

library(ggsci)
p = ggplot(plot.data.sub)+
  geom_point(aes(x = log10(pthres),y = result,colour =eth))+
  geom_line(aes(x = log10(pthres),y = result,colour = eth))+
  theme_Publication()+
  scale_color_nejm()+
  xlab("log10(P-value)")+
  ylab("AUC")+
  facet_grid(cols=vars(trait))+
  labs(colour = "Ethnic Group")

png("../../qc_plot/PT_binary_trait.png",width = 12,height = 6, units = "in",res = 300)
print(p)
dev.off()


plot.data.sub = plot.data %>% 
  filter(trait %in% c( "Heart metabolic disease burden",
                       "Height"))
p = ggplot(plot.data.sub)+
  geom_point(aes(x = log10(pthres),y = result,colour =eth))+
  geom_line(aes(x = log10(pthres),y = result,colour = eth))+
  theme_Publication()+
  scale_color_nejm()+
  xlab("log10(P-value)")+
  ylab("R2")+
  facet_grid(cols=vars(trait))+
  labs(colour = "Ethnic Group")
png("../../qc_plot/PT_con_trait.png",width = 8,height = 6, units = "in",res = 300)
print(p)
dev.off()




pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
method = "tdld"
#load all results
plot.data.list = list()
library(data.table)


r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec = c(50,100)

temp = 1
total = length(r2_vec)*length(wc_base_vec)*length(pthres)*length(pthres)
r2 = rep(0,total)
wc = rep(0,total)
ptar = rep(0,total)
peur = rep(0,total)
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
for(k1 in 1:length(pthres)){
  for(k2 in 1:length(pthres)){
    r2[temp] = r2_vec[r_ind]
    wc[temp] = wc_vec[w_ind]
    ptar[temp] = pthres[k1]
    peur[temp] = pthres[k2]
    temp = temp+1
  }
}
  }
}

index.table = data.frame(r2,wc,ptar,peur)
result.list = list()
temp = 1
for(l in 1:7){
  for(i in 2:5){
    result = read.csv(paste0(eth_group[i],"_",trait[l],"_",method))
    best.ind = index.table[which.max(result),]
    result = data.frame(eth_name = eth_name[i],
                        trait_name = trait_name[l],
                        best.ind)
    result.list[[temp]] = result
    temp = temp +1 
    # if(length(result)!=972){
    #   print(c(l,i))
    # }
    # plot.data = data.frame(pthres = pthres[(1+length(pthres)-length(result)):length(pthres)],
    #                        result = as.numeric((result)),
    #                        eth = eth_name[i],
    #                        trait = trait_name[l])
    # 
    # plot.data.list[[temp]] = plot.data
    # temp = temp+1
  }
}
result = rbindlist(result.list)
write.csv(result,file = "../../qc_plot/TDLD_bestcutoff.csv",row.names = F,quote = F)
