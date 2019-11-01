#prepare plot for haplotype view to plot
library(data.table)


#read the legend file
leg <- as.data.frame(fread("/spin1/users/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr22.legend",header=T))
leg_hap <- as.data.frame(fread("/spin1/users/zhangh24/multi_ethnic/result/LD_simulation/EUR.legend",header=T))
colnames(leg_hap)[1] <- "id"
library(dplyr)

leg_hap = left_join(leg_hap,leg,by="id")


dim(leg)
dim(leg_hap)
#merge the leg and hap file into haps file
#the output of hapgen2 called it haps, but it's actually hap file
#only take the 10,000 rows for plotting purpose
#paste -d" " EUR_noheader.legend EUR.controls.haps | head -10000 > EUR.controls_combine.haps 
haps <- as.data.frame(fread("/spin1/users/zhangh24/multi_ethnic/result/LD_simulation/EUR.controls_combine.haps",header=F))
#only take the Biallelic_SNP for plot
idx <- which(leg_hap$TYPE=="Biallelic_SNP")
idx <- idx[idx<=10000]
haps_plot = haps[idx,-c(1,2)]
fwrite(haps_plot,file = paste0("/spin1/users/zhangh24/multi_ethnic/result/LD_simulation/EUR.controls_plot.haps"),col.names = F,
       row.names=F,sep=" ",quote=F)





library(dplyr)
leg_temp = leg %>% 
  filter(AMR>=0.005&
           EAS>=0.005&
           EUR>=0.005&
           SAS>=0.005&
           AFR>=0.005)
