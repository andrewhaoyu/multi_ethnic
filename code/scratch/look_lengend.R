library(data.table)

leg <- as.data.frame(fread("/spin1/users/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr22.legend",header=T))

library(dplyr)
leg_temp = leg %>% 
  filter(AMR>=0.01&
          EAS>=0.01&
           EUR>=0.01&
           SAS>=0.01&
           AFR>=0.01&TYPE=="Biallelic_SNP") %>% 
  select(position)


idx <- which(duplicated(leg_temp))
write.table(leg_temp,file = paste0("/data/zhangh24/KG.impute2/tag/chr",22,".tag"),row.names = F,col.names = F,quote=F)


MyHeatmap <- LDheatmap(CEUSNP, CEUDist, LDmeasure="r",
                        title="Pairwise LD in r^2", add.map=TRUE,
                        SNP.name=c("rs2283092", "rs6979287"),
                        color=grey.colors(20), name="myLDgrob",
                        add.key=TRUE)
LDheatmap(CEUSNP, CEUDist, LDmeasure="r",
          title="Pairwise LD in r^2", add.map=TRUE,
          #SNP.name=c("rs2283092", "rs6979287"),
          color=grey.colors(20), 
          name="myLDgrob",
          add.key=TRUE)


data <- as.data.frame(fread("/Users/zhangh24/Desktop/1000GP_Phase3_chr22_tag_100.gen"))
snp <- t(data[,-c(1:6)])
LD <- cor(snp)^2
library(corrplot)
p1 = corrplot(LD,method="color",type="upper",diag = T,tl.cex = 1,tl.col = "white")
data2 <- as.data.frame(fread("/Users/zhangh24/Desktop/EUR.controls.tags_100.gen"))
snp2 <- t(data2[1:100,-c(1:5)])
LD2 <- cor(snp2)^2
p2 = corrplot(LD2,method="color",type="upper",diag = T,tl.cex = 1,tl.col = "white")
