#prepare the sample data for different ethnic groups for ploting
#cd /data/zhangh24/multi_ethnic/result/LD_simulation/EUR
#head -100  chr22.controls.tags.gen > chr22.controls.tags_100.gen
#similar process for all the ethnic groups

#use qctools to tranform the 1000 genome haplotype data into the .gen format
#only do this only chr 22
#need to change the hap and legend name to be the same
#cd /data/zhangh24/KG.impute2/EUR
#mv 1000GP_Phase3_chr22.legend chr22.legend

code <- rep("c",5)
eth <- c("EUR","AFR","AMR","EAS","SAS")
for(k in 1:5){
  code[k] <- paste0("/spin1/users/zhangh24/software/qctool/qctool -g /data/zhangh24/KG.impute2/",eth[k],"/chr22.hap -filetype impute_haplotypes -og /data/zhangh24/KG.impute2/",eth[k],"/1000GP_Phase3_chr22.gen -incl-positions /spin1/users/zhangh24/KG.impute2/tag/chr22.tag")
  
}

write.table(code,file = paste0("/data/zhangh24/multi_ethnic/code/LD_simulation/transform_1kg_hap_to_gen.sh"),row.names = F,
            col.names = F,
            quote=F)

#take the first 100 lines for the 1kg reference data
#cd /data/zhangh24/KG.impute2/EUR
#head -100 1000GP_Phase3_chr22.gen > 1000GP_Phase3_chr22_100.gen

#download 1kg data and simulated data to local /Users/zhangh24/GoogleDrive/multi_ethnic/result/LD_simulation/
#load data for plot
library(corrplot)
eth <- c("EUR","AFR","AMR","EAS","SAS")
setwd("/Users/zhangh24/GoogleDrive")
for(i in 1:5){
  data <- as.data.frame(fread(paste0("./multi_ethnic/result/LD_simulation/",eth[i],"/1000GP_Phase3_chr22_100.gen")))
  snp <- t(data[,-c(1:6)])
  LD <- cor(snp)^2
  png(file = paste0("./multi_ethnic/result/LD_simulation/",eth[i],"/1kg_heatmap.png"),width=6,height=6,units="in",res=300)
  corrplot(LD,method="color",type="upper",diag = T,tl.cex = 1,tl.col = "white")
  dev.off()
  data2 <- as.data.frame(fread(paste0("./multi_ethnic/result/LD_simulation/",eth[i],"/chr22.controls.tags_100.gen")))
  snp2 <- t(data2[1:100,-c(1:5)])
  LD2 <- cor(snp2)^2
  png(file = paste0("./multi_ethnic/result/LD_simulation/",eth[i],"/simulated_heatmap.png"),width=6,height=6,units="in",res=300)
  corrplot(LD2,method="color",type="upper",diag = T,tl.cex = 1,tl.col = "white")
  dev.off()

  }
