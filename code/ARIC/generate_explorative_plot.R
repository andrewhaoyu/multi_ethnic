eth <- c("EUR","AFR","AMR","EAS","SAS")
trait = c("eGFRcr","ACR","urate")

library(data.table)
y.list <- list()

i = 2
l = 1
temp = 1
for(i in 1:2){
  for(l in c(1,3)){
    
    data.dir = "/Users/zhangh24/GoogleDrive/multi_ethnic/data/ARIC/"
    
    
    #y <- as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/pheno/pheno.txt")))
    y <- as.data.frame(fread(paste0(data.dir,trait[l],"/",eth[i],"/pheno.txt")))
    colnames(y)[2] = "ID"
    eth.vec = rep(eth[i],nrow(y))
    trait.vec = rep(trait[l],nrow(y))
    y$eth.vec = eth.vec
    y$trait.vec = trait.vec
    y.list[[temp]] = y
    temp = temp+1
  }
}
y.plot = rbindlist(y.list)
#plot(y$pc1,y$pc2)
library(ggplot2)
library(dplyr)
y.plot.filter = y.plot %>% 
  filter(trait.vec == "eGFRcr") %>% 
  mutate(sex.var = ifelse(sex==1,"Male","Female")) %>% 
  mutate(sexeth = paste0(sex.var,"-",eth.vec))
ggplot(y.plot.filter,aes(pc1,pc2,col=eth.vec))+
  geom_point()+
  theme_Publication()+
  scale_colour_Publication()




ggplot(y.plot.filter,)+
  geom_density(alpha=.5,aes(x = y,y=..density..,fill=sexeth)) +
  theme_Publication()+
  scale_fill_Publication()+
  ggtitle("eGFRcr")


y.plot.filter = y.plot %>% 
  filter(trait.vec == "eGFRcr") %>% 
  mutate(sex.var = ifelse(sex==1,"Male","Female")) %>% 
  mutate(sexeth = paste0(sex.var,"-",eth.vec))

y.plot.filter = y.plot %>% 
  filter(trait.vec == "eGFRcr") %>% 
  mutate(eth.vec=="AFR")


  mutate(sex.var = ifelse(sex==1,"Male","Female")) %>% 
  mutate(sexeth = paste0(sex.var,"-",eth.vec))

ggplot(y.plot.filter)+
  geom_histogram(alpha=.5,aes(x = y)) +
  theme_Publication()+
  scale_fill_Publication()+
  ggtitle("eGFRcr")

y.plot.filter = y.plot %>% 
  filter(trait.vec == "urate")%>% 
  mutate(sex.var = ifelse(sex==1,"Male","Female")) %>% 
  mutate(sexeth = paste0(sex.var,"-",eth.vec))

ggplot(y.plot.filter)+
  geom_density(alpha=.5,aes(x = y,y=..density..,fill=sexeth)) +
  theme_Publication()+
  scale_fill_Publication()+
  ggtitle("urate")




ggplot(y.plot.filter)+
  geom_density(alpha=.5,aes(x = age,y=..density..,fill=sexeth)) +
  theme_Publication()+
  scale_fill_Publication()+
  ggtitle(".")


data = read.csv(paste0(data.dir,"/eGFRcr/AFR/egfr_pc_AA_withgwasid.csv"),header=T)

ggplot(data)+
  geom_histogram(alpha=.5,aes(x = egfrcr_v1)) +
  theme_Publication()+
  scale_fill_Publication()+
  ggtitle(".")

