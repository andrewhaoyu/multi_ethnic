source("/data/zhangh24/multi_ethnic/code/plot_code/theme_publication.R")
#load the SNPs information
load("/data/zhangh24/multi_ethnic/result/LD_simulation_new/snp.infor.rdata")
library(tidyverse)

snp.infor = snp.infor %>% 
  mutate(EUR.bi = ifelse(EUR >=0.01&EUR<=0.99,TRUE,FALSE),
         AFR.bi = ifelse(AFR >=0.01&AFR<=0.99,TRUE,FALSE),
         AMR.bi = ifelse(AMR >=0.01&AFR<=0.99,TRUE,FALSE),
         EAS.bi = ifelse(EAS >= 0.01& EAS <=0.99,TRUE,FALSE),
         SAS.bi = ifelse(SAS >= 0.01& EAS <=0.99,TRUE,FALSE))

bi.mat = snp.infor %>% select(EUR.bi,AFR.bi,
                              AMR.bi,EAS.bi,SAS.bi) 
eth = c("EUR","AFR","AMR","EAS","SAS")
total = length(eth)
eth.vec = rep("c",total)
n.specific = rep(0,total)
n.sharedEA = rep(0,total)
n.sharednonEA = rep(0,total)







for(i in 1:length(eth)){
  if(i==1){
    idx = which(bi.mat[,1]==1&
                  rowSums(bi.mat[,2:5])==0)
    n.specific[i] = length(idx)
    idx = which(bi.mat[,1]==1&
                  rowSums(bi.mat[,2:5])!=0)
    n.sharednonEA[i] = length(idx)
    n.sharedEA[i] = 0
    
  }else{
    other.ind = setdiff(c(1:5),c(i))
    nonEA.ind = setdiff(c(1:5),c(1,i))
    idx = which(bi.mat[,i]==1&
                  rowSums(bi.mat[,other.ind])==0)
    n.specific[i] = length(idx)
    idx = which(bi.mat[,1]==1&
                  bi.mat[,i]==1)
    n.sharedEA[i] = length(idx)
    
    idx = which(bi.mat[,1]==0&
                  bi.mat[,i]==1&
                  rowSums(bi.mat[,nonEA.ind])!=0)
    n.sharednonEA[i] = length(idx)
    
    
  }
}
n.specific + n.sharednonEA + n.sharedEA
colSums(bi.mat)


result.data <- data.frame(eth,n.specific,
                          n.sharedEA,
                          n.sharednonEA)
colnames(result.data) = c("Ethnic Group",
                          "Population-specific SNPs",
                          "SNPs shared with EUR",
                          "SNPs only shared with non-EUR ethnic groups")
save(result.data, file = "/data/zhangh24/multi_ethnic/result/LD_simulation_new/SNPs_distribution.rdata")
library(ggplot2)
library(reshape2)
load("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA/SNPs_distribution.rdata")
library(tidyverse)
total = rowSums(result.data[,2:4])
colnames(result.data) = c("Ethnic group",
                          "Population-specific SNPs",
                          "SNPs shared with EA",
                          "SNPs only shared with non-EA populations")
result.pro = result.data[,2:4]/total

pro.result = cbind(eth = result.data[,1],result.pro)
data_long <- melt(pro.result, id.vars=c("eth"),
                  variable.name = "SNP_Categories",
                  value.name = "Proportions") 
colnames(data_long)[1] = c("eth")
  data_long = data_long %>%
    #filter(eth!="EUR") %>%
    mutate(eth = factor(eth,levels = c("EUR","AFR","AMR","EAS","SAS")),
                                   SNP_Categories = as.character(SNP_Categories))
  data_long[6,2] = "EA SNPs shared with other non-EA populaitons"
  data_long = data_long %>% mutate(SNP_Categories = factor(SNP_Categories,
                                                              levels = c("SNPs shared with EUR",
                                                                          "SNPs only shared with non-EUR populations",
                                                                         "EA SNPs shared with other non-EUR populaitons",
                                                                          "Population-specific SNPs")))
scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#ef3b2c","#386cb0","#EF7E3D","#FFD042","#7fc97f","#004953","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

p = ggplot(data_long,aes(eth,Proportions,fill = SNP_Categories))+
  geom_bar(stat="identity")+
  theme_Publication()+
  #scale_fill_brewer(palette="Blues")+
  scale_fill_manual(values=c("#EF7E3D","#FFD042","#999999","#386cb0"))+
  ylab("Proportions")+
  xlab("Ethnic Group")+
  labs(fill='SNP Categories')+
  ggtitle(" Proportions for population-specific SNPs and shared SNPs")

png(filename = "/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA/SNPs_distribution_across_eth.png",units = "in",width = 10, height = 8,res = 300)
print(p)
dev.off()

load("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA/SNPs_distribution.rdata")
library(tidyverse)
total = rowSums(result.data[,2:4])
colnames(result.data) = c("Ethnic group",
                          "Population-specific SNPs",
                          "SNPs shared with EA",
                          "SNPs only shared with non-EA populations")
result.pro = result.data[,2:4]/total
result.pro.shared = rowSums(result.pro[,2:3])
result.pro.new = cbind(result.pro[,1],result.pro.shared)
pro.result = data.frame(eth = result.data[,1],result.pro.new)
colnames(pro.result)[2:3] = c("Population-specific SNPs","SNPs shared with other ethnic groups")
data_long <- melt(pro.result, id.vars=c("eth"),
                  variable.name = "SNP_Categories",
                  value.name = "Proportions") 
colnames(data_long)[1] = c("eth")
data_long = data_long %>%
  #filter(eth!="EUR") %>%
  mutate(eth = factor(eth,levels = c("AFR","AMR","EAS","EUR","SAS")),
         SNP_Categories = as.character(SNP_Categories))
data_long = data_long %>% 
  mutate(SNP_Categories = factor(SNP_Categories,
                              levels = c("SNPs shared with other ethnic groups","Population-specific SNPs")))
p = ggplot(data_long,aes(eth,Proportions,fill = SNP_Categories))+
  geom_bar(stat="identity")+
  theme_Publication()+
  #scale_fill_Publication()+
  scale_fill_manual(values=c("#EF7E3D","#386cb0"))+
  ylab("Proportions")+
  xlab("Ethnic Group")+
  labs(fill='SNP Categories')+
  ggtitle(" Proportions for population-specific SNPs and shared SNPs")

png(filename = "/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_GA/SNPs_distribution_across_eth_simple.png",units = "in",width = 10, height = 8,res = 300)
print(p)
dev.off()

                                 