#plot the MAF freq between EUR and AFR in KG pruned
load("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/multi_ethnic/result/pruned_MAF.Rdata")
plot(pruned.snp.clean$MAF.EUR,pruned.snp.clean$all.snp.AFR.MAF.AFR)
hist(pruned.snp.clean$MAF.EUR)
hist(pruned.snp.clean$MAF.AFR)

#colnames(pruned.snp.clean)[7] <- "MAF.AFR"
library(ggplot2)
png("./multi_ethnic/result/MAF_pruned_EUR.png",width=30,
    height= 20, units="cm", res= 300)
ggplot(pruned.snp.clean,aes(x=MAF.EUR))+
  geom_bar(stat="bin",fill="dodgerblue4")+
  theme_minimal()+
  xlab("MAF of SNPs in European population")
dev.off()
png("./multi_ethnic/result/MAF_pruned_AFR.png",width=30,
    height= 20, units="cm", res= 300)
ggplot(pruned.snp.clean,aes(x=MAF.AFR))+
  geom_bar(stat="bin",fill="#c0392b")+
  theme_minimal()+
  xlab("MAF of SNPs in African population")
dev.off()
png("./multi_ethnic/result/MAF_pruned_AMR.png",width=30,
    height= 20, units="cm", res= 300)
ggplot(pruned.snp.clean,aes(x=MAF.AMR))+
  geom_bar(stat="bin",fill="chartreuse4")+
  theme_minimal()+
  xlab("MAF of SNPs in American population")
dev.off()






load("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/multi_ethnic/result/all_snp_MAF.Rdata")
png("./multi_ethnic/result/MAF_EUR.png",width=30,
    height= 20, units="cm", res= 300)
ggplot(all.snp.update,aes(x=MAF.EUR))+
  geom_bar(stat="bin",fill="dodgerblue4")+
  theme_minimal()+
  xlab("MAF of SNPs in European population")
dev.off()
png("./multi_ethnic/result/MAF_AFR.png",width=30,
    height= 20, units="cm", res= 300)
ggplot(all.snp.update,aes(x=MAF.AFR))+
  geom_bar(stat="bin",fill="#c0392b")+
  theme_minimal()+
  xlab("MAF of SNPs in African population")
dev.off()
png("./multi_ethnic/result/MAF_AMR.png",width=30,
    height= 20, units="cm", res= 300)
ggplot(all.snp.update,aes(x=MAF.AMR))+
  geom_bar(stat="bin",fill="chartreuse4")+
  theme_minimal()+
  xlab("MAF of SNPs in American population")
dev.off()


