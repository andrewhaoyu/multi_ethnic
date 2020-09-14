#create venn diagram for different ethnic groups
load("/Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation_new/cau.snp.infor.list.rdata")
EUR<- rep(0,3)
AFR <- rep(0,3)
AMR <- rep(0,3)
EAS <- rep(0,3)
SAS <- rep(0,3)
l.vec <- rep(0,3)
for(l in 1:3){
  
  cau.snp <- cau.snp.infor.list[[l]]
  l.vec[l] <- l
  
  EUR.ind <- which(cau.snp$EUR>=0.01&
                 cau.snp$EUR<=0.99)
  EUR[l] <- length(EUR.ind)
  AFR.ind <-  which(cau.snp$AFR>=0.01&
                  cau.snp$AFR<=0.99)
  AFR[l] <- length(AFR.ind)
  AMR.ind <-  which(cau.snp$AMR>=0.01&
                  cau.snp$AMR<=0.99)
  AMR[l] <- length(AMR.ind)
  EAS.ind <-  which(cau.snp$EAS>=0.01&
                  cau.snp$EAS<=0.99)
  EAS[l] <- length(EAS.ind)
  SAS.ind <-  which(cau.snp$SAS>=0.01&
                  cau.snp$SAS<=0.99)
  SAS[l] <- length(SAS.ind)
  
}

library(tidyr)
data <- data.frame(l.vec,EUR,AFR,AMR,
              EAS,SAS)
data.long <- gather(data,ethnic,number,EUR:SAS)

data.temp <- data.long[data.long$l.vec==2,]

library(ggplot2)
data.temp$ethnic <- factor(data.temp$ethnic,
                           levels=c("EUR","AFR","AMR","EAS","SAS"))
p <- ggplot(data.temp,aes(x=ethnic,y=number))+
  geom_bar(aes(fill=ethnic),
           stat="identity")+
  ylab("Number of causal SNP")+
  xlab("Ethnic group")+
  labs(fill="Ethnic group")+
  theme_Publication()+
  scale_fill_nejm()+
  ggtitle("Number of causal SNPs (causal SNPs proportion =0.001)")
p
png(file = paste0("./Number_of_causal_SNPs_summary.png"),
    width = 10, height = 8, res = 300,units = "in")
p
dev.off()





