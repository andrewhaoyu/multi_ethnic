#create tag file for LD simulation
#tag SNPs are Biallelic_SNP with MAF 1% in EUR, AFR, AMR, EAS and SAS
#tag files are created seperately for each chr
library(data.table)
library(dplyr)
eth <- c("AFR","AMR","EAS","EUR","SAS")
for(i in 1:22){
  leg <- as.data.frame(fread(paste0("/data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",i,".legend"),header=T))
  TYPE = leg$TYPE
  MAF <- leg[,c(6:10)]
  for(k in 1:length(eth)){
    idx <- which(MAF[,k]>=0.01&
                   MAF[,k]<=0.99&
                   TYPE=="Biallelic_SNP")
    leg_chr = leg$position[idx] 
    
    
    write.table(leg_chr,file = paste0("/data/zhangh24/KG.impute2/tag/",eth[k],"_chr",i,".tag"),row.names = F,col.names = F,quote=F)
    
  }
  
  
}
