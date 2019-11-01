#create tag file for LD simulation
#tag SNPs are Biallelic_SNP with MAF 1% in EUR, AFR, AMR, EAS and SAS
#tag files are created seperately for each chr
library(data.table)
library(dplyr)
for(i in 1:22){
 
  leg <- as.data.frame(fread(paste0("/data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",i,".legend"),header=T))
  
  
  leg_chr = leg %>% 
    filter(AMR>=0.01&
             EAS>=0.01&
             EUR>=0.01&
             SAS>=0.01&
             AFR>=0.01&TYPE=="Biallelic_SNP") %>% 
    select(position)
  

  write.table(leg_chr,file = paste0("/data/zhangh24/KG.impute2/tag/chr",i,".tag"),row.names = F,col.names = F,quote=F)
  
}

